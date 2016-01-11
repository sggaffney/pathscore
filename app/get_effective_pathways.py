"""Runs pathway pipeline on CancerDB tables or TCGA tables."""

from flask import current_app
import MySQLdb as mdb
import numpy as np
from scipy import stats
from scipy.optimize import minimize_scalar, brentq
import timeit
from collections import OrderedDict
import os
import subprocess
import pyximport
pyximport.install(setup_args={'include_dirs': np.get_include()})
from sqlalchemy.orm.exc import StaleDataError
from comb_functions import get_pway_likelihood_cython

from . import db, celery
import emails
from .models import UserFile
import app
from .db_lookups import lookup_path_sizes, lookup_background_size, \
    lookup_patient_counts, lookup_patient_lengths, build_path_patient_dict, \
    lookup_path_lengths, fetch_path_ids_interest_genes, get_pathway_name_dict, \
    get_gene_combs_hit, get_gene_counts, get_pway_lenstats_dict, \
    fetch_path_info_global, count_patients, lookup_hypermutated_patients
import misc
import naming_rules

# GLOBALS
path_size_dict_full = lookup_path_sizes()
path_len_dict_full = lookup_path_lengths()
path_info_dict = fetch_path_info_global()

background_len_full = lookup_background_size(use_length=True)
background_count_full = lookup_background_size(use_length=False)


class TableLoadException(Exception):
    pass


class MatlabFailureException(Exception):
    pass


# @async
@celery.task
def run_analysis_async(proj_dir, data_path, upload_id):
    """Asynchronous run of pathway analysis."""
    user_upload = UserFile.query.get(upload_id)
    user_upload.is_queued = 0
    db.session.add(user_upload)
    db.session.commit()
    # global dbvars
    # with app.app_context():
    table = None
    unused_gene_path = naming_rules.get_unused_gene_path(user_upload)
    rejected_gene_path = naming_rules.get_rejected_gene_path(user_upload)
    try:
        # self.update_state(state='PROGRESS', meta={'status': 'Loading data'})
        table_name = user_upload.get_table_name()
        table = MutationTable(table_name, data_path,
                              unused_path=unused_gene_path,
                              rejected_path=rejected_gene_path)
        user_upload.n_rejected = table.n_rejected
        user_upload.n_ignored = table.n_ignored
        user_upload.n_loaded = table.n_loaded
        db.session.add(user_upload)
        db.session.commit()
        if not table.loaded:
            raise TableLoadException("Failed to load table {!r}."
                                     .format(table_name))
        user_upload.n_patients = count_patients(table_name)
        db.session.add(user_upload)
        db.session.commit()
        run(proj_dir, table_name, user_upload)
        user_upload.run_complete = True
    except MatlabFailureException as e:
        user_email = user_upload.uploader.email
        current_app.logger.error("Matlab failure in proj {} ({}): {}"
                                 .format(upload_id, user_email, e.args[0]))
    except TableLoadException as e:
        user_email = user_upload.uploader.email
        current_app.logger.error("Table load failure in proj {} ({}): {}"
                                 .format(upload_id, user_email, e.args[0]))
    except mdb.Error:
        pass  # already logged by db commands
    except Exception as e:  # catch all remaining errors
        current_app.logger.error("Failure in proj {} for user {}.".format(
                user_upload.file_id, user_upload.uploader.email))
    finally:
        if not user_upload.run_complete:
            user_upload.run_complete = None
            current_app.logger.info("Project {} failed for user {}.".format(
                user_upload.file_id, user_upload.uploader.email))
        else:
            current_app.logger.info("Project {} completed for user {}.".format(
                user_upload.file_id, user_upload.uploader.email))
        if table:
            drop_table(table_name)
        try:
            db.session.add(user_upload)
            db.session.commit()
            emails.run_finished_notification(upload_id)
        except StaleDataError as e:
            current_app.logger.info("Early termination of stale project: {}".format(
                e.message))


def run_analysis(proj_dir, data_path, upload_id):
    # user_upload.run_accepted = True
    # db.session.add(user_upload)
    # db.session.commit()
    # app_obj = current_app._get_current_object()
    run_analysis_async.delay(proj_dir, data_path, upload_id)


def drop_table(table_name):
    cmd = """drop table {};""".format(table_name)
    try:
        con = mdb.connect(**app.dbvars)
        cur = con.cursor()
        cur.execute(cmd)
        con.commit()
    except mdb.Error as e:
        current_app.logger.error("Error %d: %s" % (e.args[0], e.args[1]))
    finally:
        if con:
            con.close()


class NonSingleResult(Exception):
    pass


class Patient():
    def __init__(self, patient_id, n_mutated, is_mutated):
        self.patient_id = patient_id
        self.n_mutated = n_mutated
        self.is_mutated = is_mutated


class MutationTable:
    """Loads data from file into table, given table_name and file path."""
    def __init__(self, table_name, data_path, unused_path=None,
                 rejected_path=None):
        self.table_name = table_name
        self.data_path = data_path
        self.loaded = False
        self.n_initial = None
        self.n_rejected = None
        self.n_ignored = None
        self.populate_table()
        self.remove_genes_unrecognized(rejected_path)
        self.remove_genes_outwith_pathways(unused_path)
        self.n_loaded = self.n_initial - self.n_rejected - self.n_ignored
        if self.n_loaded:
            self.loaded = True

    def populate_table(self):
        """Build mutation table before filtering. Return boolean for success."""
        table_name, data_path = self.table_name, self.data_path
        create_str = u"""CREATE TABLE `{}` (
          `hugo_symbol` VARCHAR(255) DEFAULT NULL,
          `entrez_id` INT(11) DEFAULT NULL,
          `patient_id` VARCHAR(255) DEFAULT NULL,
          KEY `temp_patient_id` (`patient_id`) USING HASH,
          KEY `temp_patient_entrez` (`patient_id`,`entrez_id`) USING HASH);"""\
            .format(table_name)
        load_str = u"""load data local infile '{}'
            into table `{}` fields terminated by '\t'
            lines terminated by '\n' ignore 1 lines;""".format(
            data_path, table_name)
        cmd = u"select count(*) from `{}` m;".format(self.table_name)
        con = None
        try:
            con = mdb.connect(**app.dbvars)
            cur = con.cursor()
            cur.execute(create_str)
            cur.execute(load_str)
            cur.execute(cmd)
            result = cur.fetchone()
            self.n_initial = result[0]
            con.commit()
        except mdb.Error as e:
            current_app.logger.error("Error %d: %s" % (e.args[0], e.args[1]))
            raise
        finally:
            if con:
                con.close()

    def remove_genes_unrecognized(self, rejected_path=None):
        if not self.n_initial:
            return
        cmd1 = u"""SELECT m.* FROM `{}` m
            LEFT JOIN refs.ncbi_entrez n ON m.entrez_id = n.geneId
            WHERE n.geneId IS NULL OR m.hugo_symbol <> n.symbol;""".format(
            self.table_name)
        cmd2 = u"""delete from m using `{}` m
          LEFT JOIN refs.ncbi_entrez n ON m.entrez_id = n.geneId
            WHERE n.geneId IS NULL OR m.hugo_symbol <> n.symbol;"""\
            .format(self.table_name)
        con = None
        try:
            con = mdb.connect(**app.dbvars)
            cur = con.cursor()
            # EXPORT REJECTED GENES
            cur.execute(cmd1)
            self.n_rejected = cur.rowcount
            if self.n_rejected > 0:
                if rejected_path:
                    self.save_mutations_subset(cur, rejected_path)
                # DELETE EXTRA GENE LINES
                cur.execute(cmd2)
                con.commit()
        except mdb.Error as e:
            current_app.logger.error("Error %d: %s" % (e.args[0], e.args[1]))
            raise
        finally:
            if con:
                con.close()

    def remove_genes_outwith_pathways(self, rejected_path=None):
        """
        Remove genes outwith pathways, save to reject file,
        return remaining mutation count.

        return: remaining mutation count
        :rtype : int
        """
        cmd1 = u"""SELECT m.* FROM `{}` m
          LEFT JOIN (SELECT DISTINCT entrez_id FROM refs.pathway_gene_link) l
          ON m.entrez_id = l.entrez_id WHERE l.entrez_id IS NULL;"""\
            .format(self.table_name)
        cmd2 = u"""delete from m using `{}` m
          LEFT JOIN (SELECT DISTINCT entrez_id FROM refs.pathway_gene_link) l
          ON m.entrez_id = l.entrez_id WHERE l.entrez_id IS NULL;"""\
            .format(self.table_name)
        con = None
        try:
            con = mdb.connect(**app.dbvars)
            cur = con.cursor()
            # EXPORT EXTRA GENE LINES
            cur.execute(cmd1)
            self.n_ignored = cur.rowcount
            if self.n_ignored > 0:
                if rejected_path:
                    self.save_mutations_subset(cur, rejected_path)
                # DELETE EXTRA GENE LINES
                cur.execute(cmd2)
                con.commit()
        except mdb.Error as e:
            current_app.logger.error("Error %d: %s" % (e.args[0], e.args[1]))
            raise
        finally:
            if con:
                con.close()

    @staticmethod
    def save_mutations_subset(cursor, save_path):
        """Export mutation subset."""
        with open(save_path, 'w') as out:
            for row in cursor:
                out.write('\t'.join([str(i) for i in row]) + '\n')


class GeneMatrix:
    """Holds patient-gene matrix info for a single pathway.
    Write matrix to file with call to export_matrix."""

    def __init__(self):
        self.genePatientDict = dict()

    def add_gene_patients(self, genepatients_dict):
        """Add gene-patientList pairs into dictionary."""
        self.genePatientDict = genepatients_dict

    def export_matrix(self, outfile, exclusive_genes):
        """Writes tab-separated matrix file for patient/gene
        pairs in pathway."""
        # sort genes alphabetically then by exclusive-status then patient counts
        genesOrdered = [gene for gene in self.genePatientDict]
        genesOrdered.sort()  # sort alphabetically
        genesOrdered.sort(key=lambda gene: gene in exclusive_genes,
                          reverse=True)
        genesOrdered.sort(key=lambda gene: len(self.genePatientDict[gene]),
                          reverse=True)
        # get set of patients and count genes hit for each patient
        patient_set = set.union(
            *[set(i) for i in self.genePatientDict.values()])
        patient_counts = dict(
            zip(list(patient_set), [0] * len(patient_set)))
        for patient in patient_set:
            for gene in genesOrdered:
                if patient in self.genePatientDict[gene]:
                    patient_counts[patient] += 1
        patient_list = list(patient_set)
        patient_list.sort()  # sort alphabetically        
        patient_list.sort(key=lambda patient: patient_counts[patient],
                          reverse=True)
        outfile.write("\t".join(['GENE'] + patient_list) + '\n')
        for gene in genesOrdered:
            if gene in exclusive_genes:
                outfile.write('*')
            outfile.write(gene)
            for patient in patient_list:
                if patient in self.genePatientDict[gene]:
                    outfile.write('\t1')
                else:
                    outfile.write('\t0')
            outfile.write('\n')


class PathwaySummaryBasic():
    """Holds pathway information."""
    def __init__(self, pathway_number):
        self.path_id = pathway_number
        self.n_actual = None
        # self.patients = list() # tuples. (patient_id, n_mutations, is_mutated)
        self.n_effective = None
        self.p_value = None


class PathwaySummary(PathwaySummaryBasic):
    """Holds pathway information, and can fetch info from db."""

    def __init__(self, pathway_number, proj_abbrvs,
                 patient_ids=list(), expressed_table=None,
                 ignore_genes=list()):
        PathwaySummaryBasic.__init__(self, pathway_number)
        # above gives path_id, n_actual, n_effective, p_value
        self.patients = list()  # tuples. (patient_id, n_mutations, is_mutated)
        self.filter_patient_ids = patient_ids
        self.filter_expressed = expressed_table
        self.ignore_genes = ignore_genes
        self.proj_abbrvs = proj_abbrvs
        self.gene_coverage = OrderedDict()
        self.exclusive_genes = list()
        self.cooccurring_genes = list()
        self.geneMatrix = None
        self.runtime = None  # only set up by file reader
        self.ll_actual = None
        self.ll_effective = None
        self.D = None
        self.ne_low = None
        self.ne_high = None

    def set_up_from_file(self, pval, psize, peffect, ll_actual, ll_effective, d,
                         ne_low, ne_high, runtime):
        """Manually specify pathway summary properties (after running init)."""
        self.p_value = pval
        self.n_actual = psize
        self.n_effective = peffect
        self.ll_actual = ll_actual
        self.ll_effective = ll_effective
        self.D = d
        self.ne_low = ne_low
        self.ne_high = ne_high
        self.runtime = runtime
        # # fetch gene_coverage, exclusive_genes, cooccurring genes
        # self._populate_exclusive_cooccurring()
        # self._update_gene_coverage()

    def set_pathway_size(self, path_size_dict):
        self.n_actual = path_size_dict[self.path_id]

    def _build_patient_filter_str(self, form="WHERE"):
        """build SQL substring to filter patients by ids in mutation lookup."""
        if form == 'WHERE':
            prepend = "WHERE "
        elif form == 'AND':
            prepend = "AND "
        else:
            raise Exception("Unrecognized form in patient filter.")
        if self.filter_patient_ids:
            filter_str = (prepend + "patient_id IN " +
                      str(self.filter_patient_ids).replace("[", "(").replace(
                          "]", ")"))
        else:
            filter_str = ""
        return filter_str

    def _build_ignore_gene_filter_str(self, form="WHERE"):
        """build SQL substring to filter genes in mutation lookup.
        e.g. (WHERE/AND) hugo_symbol in ('BRAF')"""
        if form == 'WHERE':
            prepend = "WHERE "
        elif form == 'AND':
            prepend = "AND "
        else:
            raise Exception("Unrecognized form in 'ignore gene' filter.")
        if self.ignore_genes:
            filter = (prepend + "NOT hugo_symbol IN " +
                      str(self.ignore_genes).replace("[", "(").replace("]",
                                                                       ")"))
        else:
            filter = ""
        return filter

    def _build_expressed_filter_str(self, join_ref):
        """build SQL substring to filter genes by entrez_id in mutation lookup.
        join_ref is abbreviation of table and column to join to expression
        table, e.g. m.entrez_id or pwg.entrez_gene_id.
        Assumes expression table is in tcga database.
        """
        if self.filter_expressed:
            filter_str = "INNER JOIN {filter_expressed} e " \
                "ON {join_ref} = e.entrez_id" \
                .format(filter_expressed=self.filter_expressed,
                        join_ref=join_ref)
        else:
            filter_str = ""
        return filter_str

    def update_exclusive_cooccurring_coverage(self, genelists, n_patients,
                                              gene_patients):
        """POSTPROCESSING.
        Gather pathway gene info and write detailed output."""
        self.geneMatrix = GeneMatrix()  # populated during update_gene_coverage
        self._populate_exclusive_cooccurring(genelists)
        self._update_gene_coverage(n_patients, gene_patients)

    def _populate_exclusive_cooccurring(self, genelists):
        """ Postprocessing step. Look up gene combinations hit
        (via _get_gene_combs_hit_yale and _get_gene_combs_hit_tcga),
        sort into sets."""
        all_hit_genes_set = set()
        exclusive_gene_set = set()
        # get minimum set of gene combinations
        temp_set = set()
        for genes in genelists:
            temp_set.add(tuple(genes))
        gene_combs_list = [list(genes) for genes in temp_set]
        for geneList in gene_combs_list:
            if len(geneList) == 1:
                exclusive_gene_set.add(geneList[0])
            for gene in geneList:
                all_hit_genes_set.add(gene)
        cooccurring_gene_set = all_hit_genes_set.difference(exclusive_gene_set)
        self.exclusive_genes = sorted(list(exclusive_gene_set))
        # self.cooccurring_genes = sorted(tagAlongGenes)
        self.cooccurring_genes = sorted(list(cooccurring_gene_set))
        return

    def _update_gene_coverage(self, n_patients, gene_patients):
        """Populate object's gene_coverage dictionary."""
        gene_coverage = dict()
        for gene in gene_patients:
            gene_coverage[gene] = len(gene_patients[gene])

        coverage_tuple = [tuple([gene, gene_coverage[gene]]) for gene in
                          gene_coverage]
        coverage_tuple.sort(key=lambda x: x[0])  # alphabeticize
        coverage_tuple.sort(key=lambda x: x[1],
                            reverse=True)  # sort by n_patients
        gene_coverage = OrderedDict(coverage_tuple)

        # convert to coverage by dividing by number of patients
        for gene in gene_coverage:
            gene_coverage[gene] = float(gene_coverage[gene]) / n_patients * 100
        self.gene_coverage = gene_coverage
        return


class LCalculator():
    """Calculates likelihood of observing pathway mutations in patients, and 
    MLE pathway size from these observations. Takes a pathwaySummary object."""

    def __init__(self, pway, genome_size=18852):
        self.G = genome_size  # genes in genome
        self.pway = pway
        self.likelihood = None
        self.ne = None
        self.ne_ll = None
        self.D = None
        self.pvalue = None
        self.ne_low = None
        self.ne_high = None
        # self.max_mutations = pway.max_mutations

        self.n_patients, self.n_mutated_array, self.is_mutated_array = \
            self._get_patient_arrays()

    def run(self):
        """ Calculate likelihood and maximum likelihood estimate."""
        self.likelihood = self._get_pway_likelihood(self.pway.n_actual)
        (ne, lastll) = self._get_ne()
        self.ne = ne
        self.ne_ll = lastll
        upper_CI, lower_CI = self._get_ne_CI()
        self.ne_low = int(upper_CI)
        self.ne_high = int(lower_CI)
        self.D = -2 * self.likelihood + 2 * self.ne_ll
        self.pvalue = 1 - stats.chi2.cdf(self.D, 1)
        self.pway.n_effective = self.ne
        self.pway.p_value = self.pvalue

    def _get_patient_arrays(self):
        """converts list of Patient objects (with n_mutated and is_mutated
        attributes) to is_mutated boolean array and n_mutated int array.
        Returns n_patients, n_mutated, is_mutated."""
        patients = self.pway.patients
        n_patients = len(patients)
        n_mutated_array = np.array([p.n_mutated for p in patients], dtype=int)
        is_mutated_array = np.array([p.is_mutated for p in patients],
                                    dtype=int)
        return n_patients, n_mutated_array, is_mutated_array

    def _get_pway_likelihood(self, pway_size=None):
        """Calculate pathway log likelihood at stated size."""
        return get_pway_likelihood_cython(self.G, pway_size,
                                          self.n_patients, self.n_mutated_array,
                                          self.is_mutated_array)

    def _get_pway_likelihood_neg(self, pway_size=None):
        """Calculate pathway log likelihood at stated size."""
        ll = get_pway_likelihood_cython(self.G, pway_size,
                                         self.n_patients, self.n_mutated_array,
                                         self.is_mutated_array)
        ll = np.nan if ll == -np.inf else ll
        return -1 * ll

    def _get_ne_CI(self):
        # if everyone has mutation in the pathway, set CI_low and high to be G
        CI_low = None
        CI_high = None
        if 0 not in self.is_mutated_array:  # everyone has mutation in pathway
            CI_high = self.G
        if 1 not in self.is_mutated_array:  # no one has a mutation in pathway
            CI_low = 1
        if self.ne == self.G:
            CI_high = self.G

        # we know there is at least one patient without a mutation, so ne
        #     should exist, and be 1 or above.
        if CI_low is None:
            try:
                CI_low = np.floor(brentq(self._get_pway_likelihood_CI, 1, self.ne))
            except ValueError:  # if signs don't change, 1 is within CI
                CI_low = 1

        if CI_high is None:  # ne is below genome size, so find maximum
            try:
                CI_high = np.ceil(brentq(self._get_pway_likelihood_CI, self.ne, self.G - self.pway.n_actual))
            except ValueError: # if signs don't change, G is within CI
                CI_high = self.G
        return CI_low, CI_high

    def _get_pway_likelihood_CI(self, pway_size=None):
        """Calculate pathway log likelihood at stated size."""
        ll = get_pway_likelihood_cython(self.G, pway_size,
                                        self.n_patients, self.n_mutated_array,
                                        self.is_mutated_array)
        return ll - self.ne_ll + 1.92

    def _get_ne_old(self):
        last_ll = None
        # improved = False
        ne = None
        # profile = list()
        # if pathway_size is zero, effective size is zero.
        if not self.pway.n_actual:
            ne = 0
            ult = np.float64(0)
            current_app.logger.debug(
                "Pathway {} contains zero genes. ".format(self.pway.path_id))
            return ne, ult
        # if all patients mutated, use ne=genome_size - max_mutations
        if False not in [patient.is_mutated for patient in self.pway.patients]:
            ne = self.G
            ult = np.float64(0)
            # current_app.logger.debug("All patients have mutation in pathway {}. ".format(
            #     self.pway.path_id) + "Effective size is full genome.")
            return (ne, ult)
        # check last 2 vals to check for decline:
        penult = self._get_pway_likelihood(
            pway_size=self.G - self.pway.n_actual - 2)  # WAS self.G - self.max_mutations - 1
        ult = self._get_pway_likelihood(
            pway_size=self.G - self.pway.n_actual - 1)  # WAS self.G - self.max_mutations
        if ult > penult:
            ne = self.G
            return (ne, ult)
        # at this stage, there will be a max before Genome size
        for pway_size in xrange(1, self.G):
            # WAS xrange(1,self.G - self.max_mutations):
            this_ll = self._get_pway_likelihood(pway_size=pway_size)
            # profile.append(this_ll)
            # if mod(pway_size,100)==0:
            # print this_ll
            if last_ll is None:
                last_ll = this_ll
            # if improved is False and (this_ll > self.likelihood or
            # pway_size >= self.pway_size):
            # improved = True
            if this_ll < last_ll or this_ll == 0:
                ne = pway_size - 1
                if this_ll == 0:
                    current_app.logger.debug("Premature stop for pway {}.".format(
                        self.pway.path_id))
                break
            last_ll = this_ll
        return ne, last_ll

    def _get_ne(self):
        """Find maximum likelihood pathway size and corresponding
        log likelihood."""

        # INITIAL BOUNDARY TESTS
        # if pathway_size is zero, effective size is zero.
        if not self.pway.n_actual:
            ne = 0
            final_ll = np.float64(0)
            current_app.logger.debug("Pathway {} contains zero genes. ".format(
                self.pway.path_id))
            return ne, final_ll
        # if all patients mutated, use ne=genome_size - max_mutations
        if False not in [patient.is_mutated for patient in self.pway.patients]:
            ne = self.G
            final_ll = np.float64(0)
            current_app.logger.debug("All patients have mutation in pathway {}. ".format(
                self.pway.path_id) + "Effective size is full genome.")
            return ne, final_ll

        result = minimize_scalar(self._get_pway_likelihood_neg, method='Brent',
                                     bounds=[1, self.G])
        # res2 = minimize_scalar(self._get_pway_likelihood_neg, method='Bounded',
        #                        bounds=[1, 17000000])
        try_integers = [np.floor(result.x), np.ceil(result.x)]
        out = [self._get_pway_likelihood_neg(i) for i in try_integers]
        if out[1] < out[0]:
            use_ind = 1
        else:
            use_ind = 0
        ne = int(try_integers[use_ind])
        final_ll = -1 * out[use_ind]
        return ne, final_ll


class GenericPathwayFileProcessor():
    """Generic object that can convert yale_proj_ids and tcga_proj_abbrvs
    to file_name."""
    base_str = 'pathways_pvalues'

    def __init__(self, dir_path, file_id, name_suffix=None):
        self.file_id = file_id
        self.dir_path = dir_path
        self.name_suffix = name_suffix
        self.root_name = self._get_root_filename()

    def _get_root_filename(self):
        """Root file name for output files."""
        root_name = os.path.join(self.dir_path, self.base_str)
        if self.name_suffix:
            root_name += '_' + self.name_suffix
        else:
            root_name += '_' + str(self.file_id)
        return root_name


class PathwayBasicFileWriter(GenericPathwayFileProcessor):
    """Writes initial p-value file."""

    def write_pvalue_file(self, lcalc, runtime):
        """Write initial processing file with p-value and MLE estimate."""
        bufsize = 1  # line buffered output
        outfile_name = self.root_name + '.txt'
        path_id = lcalc.pway.path_id
        with open(outfile_name, 'a', bufsize) as out:
            out.write('{}\t{:.3e}\t{}\t{}\t{:.3g}\t{:.3g}\t{:.3g}\t{}\t{}\t{:.2f}\n'.format(
                path_id, lcalc.pvalue, lcalc.pway.n_actual, lcalc.ne,
                lcalc.likelihood, lcalc.ne_ll, lcalc.D,
                lcalc.ne_low, lcalc.ne_high,
                runtime))


class PathwayListAssembler(GenericPathwayFileProcessor):
    """Builds ordered list of pathways from basic p-value file."""

    def __init__(self, dir_path, file_id, proj_abbrvs, patient_ids=list(),
                 name_suffix=None, expressed_table=None,
                 ignore_genes=list()):
        # create self.root_name
        GenericPathwayFileProcessor.__init__(self, dir_path, file_id,
                                             name_suffix=name_suffix)
        self.proj_abbrvs = proj_abbrvs
        self.filter_patient_ids = patient_ids
        self.expressed_table = expressed_table
        self.ignore_genes = ignore_genes

    def get_ordered_pway_list(self):
        file_name = self.root_name + '.txt'
        # max_lookup_rows = 100
        all_pathways = list()
        with open(file_name, 'r') as file_in:
            for line in file_in:
                row = line.strip().split()
                path_id = int(row[0])
                pval = float(row[1])
                psize = int(row[2])
                peffect = int(row[3])
                ll_actual = float(row[4])
                ll_effective = float(row[5])
                D = float(row[6])
                ne_low = int(row[7])
                ne_high = int(row[8])
                runtime = float(row[9])
                # set up pathway object
                pway = PathwaySummary(path_id, self.proj_abbrvs,
                                      patient_ids=self.filter_patient_ids,
                                      expressed_table=self.expressed_table,
                                      ignore_genes=self.ignore_genes)
                pway.set_up_from_file(pval, psize, peffect, ll_actual,
                                      ll_effective, D, ne_low, ne_high, runtime)
                all_pathways.append(pway)
        # sort pathways by p-value FIRST
        all_pathways.sort(key=lambda pw: pw.D, reverse=True)
        # NEXT sort pathways by effect_size_lowCI : actual_size
        all_pathways.sort(key=lambda p:
                          self._get_size_ratio(p.n_effective, p.n_actual), reverse=True)

        # # for pway in allPathways[0:max_lookup_rows+1]:
        # for pway in allPathways:
        # if pway.p_value < 0.1:
        # pway.update_exclusive_cooccurring_coverage()
        return all_pathways

    @staticmethod
    def _get_size_ratio(n_effective, n_actual):
        """Ratio used for initial sorting. >1 if large effective.
        0<r<1 if small_effective. 0 if n_actual is 0."""
        if n_actual:
            return float(n_effective) / n_actual
        else:
            return 0


class PathwayDetailedFileWriter(GenericPathwayFileProcessor):
    """Writes detailed postprocessing file:
    pathway names, pvalues and gene info."""

    name_postfix = '_detail.txt'
    def __init__(self, dir_path, file_id, pway_object_list,
                 name_suffix=None, path_genelists_dict=dict(),
                 path_genepatients_dict=dict(), n_patients=None):
        # create self.root_name
        GenericPathwayFileProcessor.__init__(self, dir_path, file_id,
                                             name_suffix=name_suffix)
        self.allPathways = pway_object_list
        self.nameDict = get_pathway_name_dict()
        self.outfile_name = self.root_name + self.name_postfix
        self.matrix_folder = 'matrix_txt'
        self.path_genelists_dict = path_genelists_dict
        self.path_genepatients_dict = path_genepatients_dict
        self.n_patients = n_patients

    @staticmethod
    def dict_to_struct(coverage_dict):
        """ Get matlab command to convert dictionary to structure.
        e.g. Dict: 'BRAF' -> 34. --> struct('BRAF',34)."""
        # struct_string = repr(coverage_dict) # "('BRAF':34, 'KRAS':16)"
        # struct_string = struct_string.replace(":",",")
        # struct_string = struct_string.replace("{","(").replace("}",")")
        # struct_string = "struct" + struct_string
        # return struct_string
        # name,coverage pairs
        pairs = ",".join(["{!r},{:.2f}".format(gene, coverage_dict[gene]) for
                          gene in coverage_dict])
        struct_string = "struct(" + pairs + ")"
        return struct_string

    def write_detailed_file(self):
        """Perform lookup of coverage etc for low pvalue pathways and write 
        ordered list of pathways plus info to file."""
        bufsize = 1
        with open(self.outfile_name, 'w', bufsize) as out:
            for pway in self.allPathways:
                # Get extra info, if p_value is low
                if pway.p_value < 0.05 and pway.n_actual < pway.n_effective \
                        and pway.n_actual < pway.ne_low:
                    pway.update_exclusive_cooccurring_coverage(
                        self.path_genelists_dict[pway.path_id],
                        self.n_patients,
                        self.path_genepatients_dict[pway.path_id])
                    if pway.gene_coverage:  # if genes are hit...
                        if pway.n_effective > pway.n_actual:
                            pway.geneMatrix.add_gene_patients(
                                self.path_genepatients_dict[pway.path_id])
                            self.write_matrix_files(pway)
                path_name = self.nameDict[pway.path_id]
                coverage_string = self.dict_to_struct(pway.gene_coverage)
                pway.exclusive_genes.sort(key=lambda gene:
                                          pway.gene_coverage[gene],
                                          reverse=True)
                pway.cooccurring_genes.sort(key=lambda gene:
                                            pway.gene_coverage[gene],
                                            reverse=True)
                exclusive_string = '{' + ','.join(
                    [repr(i) for i in pway.exclusive_genes]) + '}'
                cooccurring_string = '{' + ','.join(
                    [repr(i) for i in pway.cooccurring_genes]) + '}'
                format_str = "{path_id}\t{name}\t{n_actual}\t{n_effective}\t" \
                             "{p_value:.3e}\t{ll_actual:.4g}\t{ll_effective:.4g}\t" \
                             "{D:.4g}\t{ne_low}\t{ne_high}\t{runtime:.2f}\t{exclusive_string}\t" \
                             "{cooccurring}\t{coverage_string}\t{pway_ngenes}\n"""
                out.write(format_str.format(path_id=pway.path_id,
                                            name=path_name,
                                            coverage_string=coverage_string,
                                            n_actual=pway.n_actual,
                                            n_effective=pway.n_effective,
                                            p_value=pway.p_value,
                                            ll_actual=pway.ll_actual,
                                            ll_effective=pway.ll_effective,
                                            D=pway.D,
                                            ne_low=pway.ne_low,
                                            ne_high=pway.ne_high,
                                            runtime=pway.runtime,
                                            exclusive_string=exclusive_string,
                                            cooccurring=cooccurring_string,
                                            pway_ngenes=path_size_dict_full[pway.path_id]))

    def write_matrix_files(self, pway):
        """Write text file containing presence matrix for patient-gene pair."""
        matrix_path = os.path.join(self.dir_path, self.matrix_folder)
        if not os.path.exists(matrix_path):
            os.mkdir(matrix_path)
        matrix_filename = os.path.join(matrix_path, 'matrix_' + str(pway.path_id) + '.txt')
        with open(matrix_filename, 'w') as outfile:
            pway.geneMatrix.export_matrix(outfile, pway.exclusive_genes)


class PathwaySummaryParsed(PathwaySummaryBasic):
    """Holds pathway information, pulled from anaysis output."""
    def __init__(self, pathway_number):
        PathwaySummaryBasic.__init__(self, pathway_number)
        # above adds path_id, n_actual, n_effective, p_value
        pathway_number = int(pathway_number)  # ensure int
        self.name = None
        self.root_url = "http://www.broadinstitute.org/gsea/msigdb/cards/"
        self.url = path_info_dict[pathway_number]['url'].split(self.root_url)[1]
        self.description = path_info_dict[pathway_number]['desc']
        self.contrib = path_info_dict[pathway_number]['contrib']
        self.gene_set = set()
        self.gene_pc = dict()
        self.lengths_tuple = tuple()  # set externally

    @property
    def nice_name(self):
        name_str = misc.strip_contributors(self.name)
        name_str = name_str.replace('_', ' ')
        name_str = name_str.replace('TEL PATHWAY', 'TELOMERASE PATHWAY')
        name_str = name_str.replace('RNA PATHWAY', 'PKR SIGNALING PATHWAY')
        return name_str

    def as_string_js(self):
        """Return string for javascript."""
        outstr = "{ id:" + self.path_id + ", " + "name:'" + misc.html_quotes(self.nice_name) \
                 + "', pval:'" + self.p_value + "', size:" + str(self.n_actual) \
                 + ", effective:" + str(self.n_effective) + ", url:'" \
                 + self.url + "', contrib: '" + misc.html_quotes(self.contrib) + "' , brief: '" \
                 + misc.html_quotes(self.description) + "', lengths:" + str(list(self.lengths_tuple)) \
                 + ", geneSet: "
        if self.gene_set:
            gene_set_str = "','".join(self.gene_set)
            outstr = outstr + "['" + gene_set_str + "']}"
        else:
            outstr = outstr + "[]}"
        return outstr

    @property
    def full_url(self):
        return self.root_url + self.url

    def as_string_readable(self):
        """Return string for javascript."""
        name_str = self.name.replace('_', ' ')
        name_str = name_str.replace('BIOCARTA ', '')
        name_str = name_str.replace('REACTOME ', '')
        name_str = name_str.replace('KEGG ', '')
        name_str = name_str.replace('TEL PATHWAY', 'TELOMERASE PATHWAY')
        name_str = name_str.replace('RNA PATHWAY', 'PKR SIGNALING PATHWAY')
        outstr = "{ id:" + self.path_id + ", " + "name:'" + name_str \
                 + "', pval:'" + self.p_value + "', size:" + str(self.n_actual) \
                 + ", effective:" + str(self.n_effective) + ", url:'" \
                 + self.url + "', brief:'" + self.description + "', contrib:'" \
                 + self.contrib + "', geneSet: "
        if self.gene_set:
            gene_set_str = "','".join(self.gene_set)
            outstr = outstr + "['" + gene_set_str + "']}"
        else:
            outstr = outstr + "[]}"
        return outstr


def get_patient_list(path_id, patient_size_dict, path_patient_dict):
    """Get list of patient objects for given pathway using dictionary
    path_patient_dict."""
    patient_list = list()
    if path_id not in path_patient_dict:
        raise Exception("Requested path_id that is not in dictionary.")
    patient_set = path_patient_dict.pop(path_id)
    for patient_id in patient_size_dict:
        n_mutated = patient_size_dict[patient_id]
        is_mutated = patient_id in patient_set
        patient = Patient(patient_id, n_mutated, is_mutated)
        patient_list.append(patient)
    return patient_list


def run(dir_path, table_name, user_upload):
    """ EXAMPLE ARGUMENTS
    dir_path = /Users/sgg/Downloads/uploads/1/41/
    table_name = mutations_41
    alg='gene_count' --OR-- 'gene_length'
    user_upload is upload object (file_id, user_id, filename, ...etc)
    """
    alg = user_upload.algorithm  # 'gene_count' or 'gene_length'
    use_length = True if alg == 'gene_length' else False

    ignore = user_upload.ignore_genes
    ignore_list = str(ignore).split(',') if ignore else []

    genome_size = lookup_background_size(ignore_genes=ignore_list,
                                         use_length=use_length)
    # GET BACKGROUND SIZE (varies with ignore_list)
    if alg == 'gene_count':
        if ignore_list:
            genome_size = lookup_background_size(ignore_genes=ignore_list,
                                                 use_length=False)
        else:
            genome_size = background_count_full
    elif alg == 'gene_length':
        if ignore_list:
            genome_size = lookup_background_size(ignore_genes=ignore_list,
                                                 use_length=True)
        else:
            genome_size = background_len_full
    else:
        raise Exception("Invalid algorithm selection ({})".format(alg))
    # self.update_state(state='PROGRESS', meta={'status': 'Calculating L'})
    detail_path, matrix_folder = generate_initial_text_output(
        dir_path, table_name, user_upload, genome_size, alg=alg)
    # matrix_folder is typically 'matrix_txt'
    # self.update_state(state='PROGRESS', meta={'status': 'Plotting'})
    generate_plot_files(user_upload, detail_path, table_name, dir_path,
                        matrix_folder, genome_size)


def generate_initial_text_output(dir_path, table_name, user_upload, genome_size, alg='gene_count'):
    # max_mutations
    file_id = user_upload.file_id
    table_list = [table_name]  # code can iterate through list of tables
    # max_mutations = user_upload.n_cutoff or 500  # default max is 500
    ignore = user_upload.ignore_genes
    ignore_genes = str(ignore).split(',') if ignore else []

    proj_suffix = user_upload.get_local_filename()
    # # get patient list. maybe empty list.
    patient_list = []
    # if args.patients_file:
    #     for line in args.patients_file:
    #         temp_line = line.strip('\n')
    #         if temp_line.isdigit():
    #             patient_list.append(int(temp_line))  # integer if possible
    #         else:
    #             patient_list.append(temp_line)  # strings for tcga projects
    #     print("Loaded {} patients.".format(len(patient_list)))
    # else:
    #     print('No patients file provided. Using all patients.')

    # get genes of interest, if any
    if user_upload.required_genes:
        interest_genes = tuple(str(user_upload.required_genes).split(','))
    else:
        interest_genes = tuple()
    # print("Interest genes: {}".format(str(interest_genes)))

    # INCLUDE_GENES GIVES PATHWAY REQUIREMENT
    # EXCLUDE_GENES MEANS DON'T BASE PATHWAY-PATIENT PAIRS ON THIS GENE

    # LOAD ALL PATHWAY-PATIENT PAIRS INTO DICTIONARY.
    # e.g.  {1L: set(['yucrena', 'yuhoin']), ...}
    # used for building Patient list for pathway (popped down to nothing) in
    #   tandem with
    path_patient_dict = build_path_patient_dict(table_name, ignore_genes)
    # CREATE LIST OF PATHWAYS TO INSPECT
    #  - STRIPPED DOWN TO THOSE WITH INTEREST GENES, IF SPECIFIED
    all_path_ids = {id for id in path_patient_dict}  # all altered pathways set
    if interest_genes:
        interest_gene_path_ids = fetch_path_ids_interest_genes(interest_genes)
        all_path_ids = all_path_ids.intersection(set(interest_gene_path_ids))
        # delete unnecessary path_ids if there are interest_genes
        unwanted_ids = list()
        for path in path_patient_dict:
            if path not in all_path_ids:
                unwanted_ids.append(path)
        for i in unwanted_ids:
            path_patient_dict.pop(i)
    # create sorted list of path_ids
    all_path_ids = sorted(list(all_path_ids))

    if alg == 'gene_count':
        if ignore_genes:
            path_size_dict = lookup_path_sizes(ignore_genes)
        else:
            path_size_dict = path_size_dict_full
    elif alg == 'gene_length':
        if ignore_genes:
            path_size_dict = lookup_path_lengths(ignore_genes)
        else:
            path_size_dict = path_len_dict_full

    if alg == 'gene_count':
        patient_size_dict = lookup_patient_counts(table_name, ignore_genes)
    elif alg == 'gene_length':
        patient_size_dict = lookup_patient_lengths(table_name, ignore_genes)

    hypermutated = lookup_hypermutated_patients(table_name)
    hyper_path = naming_rules.get_hypermutated_path(user_upload)
    with open(hyper_path, 'w') as out:
        for patient_id in hypermutated:
            out.write(patient_id + '\n')

    for pathway_number in all_path_ids:
        # Populate pathway object, and time pvalue calculation
        start = timeit.default_timer()
        pway = PathwaySummary(pathway_number, table_list,
                              # max_mutations=max_mutations,
                              expressed_table=None,
                              ignore_genes=ignore_genes)
        pway.set_pathway_size(path_size_dict)
        pway.patients = get_patient_list(pathway_number, patient_size_dict,
                                         path_patient_dict)
        lcalc = LCalculator(pway, genome_size)  # include optional genome_size
        lcalc.run()
        runtime = timeit.default_timer() - start
        # Write results to 'basic' file
        basic_writer = PathwayBasicFileWriter(dir_path, file_id,
                                              name_suffix=proj_suffix)
        basic_writer.write_pvalue_file(lcalc, runtime)

    # Gather all pathway stats from text file
    assembler = PathwayListAssembler(dir_path, file_id, table_list,
                                     # patient_ids=None,
                                     name_suffix=proj_suffix,
                                     # expressed_table=args.expression,
                                     ignore_genes=ignore_genes)
    pway_list = assembler.get_ordered_pway_list()

    path_genelists_dict = get_gene_combs_hit(table_name)
    path_genepatients_dict = get_gene_counts(table_name)
    n_patients = len(patient_size_dict)

    # Rank pathways, gather extra stats and write to final file
    final_writer = PathwayDetailedFileWriter(dir_path, file_id, pway_list,
                                             name_suffix=proj_suffix,
                                             path_genelists_dict=path_genelists_dict,
                                             path_genepatients_dict=path_genepatients_dict,
                                             n_patients=n_patients)
    final_writer.write_detailed_file()

    # generate_files(dir_path, final_writer.outfile_name, user_upload)

    detail_path = final_writer.outfile_name
    # descriptive_name = user_upload.get_local_filename()
    matrix_folder = final_writer.matrix_folder
    return detail_path, matrix_folder


def generate_plot_files(user_upload, detail_path, table_name, dir_path,
                        matrix_folder, genome_size):
    allPathways = load_pathway_list_from_file(detail_path)

    scores_path, names_path, tree_svg_path = naming_rules.\
        get_tree_score_paths(user_upload)
    save_dissimilarity_files(allPathways, scores_path, names_path, min_len=50)

    # LOOKUP MUTATED GENE LENGTHS
    if user_upload.ignore_genes:
        ignore_genes = str(user_upload.ignore_genes).split(',')
    else:
        ignore_genes = []
    pathway_lengths = get_pway_lenstats_dict(table_name, ignore_genes)
    for p in allPathways:
        p.lengths_tuple = pathway_lengths[int(p.path_id)]

    js_out_path = naming_rules.get_js_path(user_upload)
    make_js_file(allPathways, js_out_path)

    readable_path = os.path.join(dir_path,
                                 user_upload.get_local_filename() + ".xls")
    make_readable_file(allPathways, readable_path)
    # html_name = create_html(detail_path, out_dir, project_str, descriptive_name,
    #                         skipfew)

    # ONLY CREATE SVGS IF MATRIX TXT PATH EXISTS (i.e. pathways have mutations)
    if os.path.exists(os.path.join(dir_path, matrix_folder)):
        hyper_path = naming_rules.get_hypermutated_path(user_upload)
        create_pway_plots(str(detail_path), scores_path, names_path,
                          tree_svg_path, genome_size, hyper_path)
        misc.zip_svgs(dir_path)
    # create_svgs(str(detail_path))
    # create_matrix_svgs(str(detail_path))


def load_pathway_list_from_file(results_path):
    """Loads detailed file into list of PathwaySummaryParsed objects.
    Used for making js file and user readable output."""
    allPathways = list()  # will hold pathway objects
    ignoreList = ['CANCER', 'GLIOMA', 'MELANOMA', 'LEUKEMIA', 'CARCINOMA']
    # project_id = raw_input("Project id? ")  # <TODO:sgg> change to input in python3
    # 'pathways_pvalues_{}_pretty.txt'.format(project_id)
    with open(results_path, 'r') as file:
        for line in file:
            vals = line.strip("\n").split('\t')
            pway = PathwaySummaryParsed(vals[0])
            # skip line if part of name is in ignoreList
            skipLine = False
            name = vals[1]
            for ignore in ignoreList:
                if ignore in name:
                    skipLine = True
                    break
            if skipLine:
                continue
            pway.name = name
            pway.n_actual = int(vals[2])
            pway.n_effective = int(vals[3])
            pway.p_value = vals[4]
            pway.ll_actual = float(vals[5])
            pway.ll_effective = float(vals[6])
            pway.D = float(vals[7])
            pway.ne_low = int(vals[8])
            pway.ne_high = int(vals[9])
            gene_set = set.union(set(eval(vals[11])), set(eval(vals[12])))
            pway.gene_set = gene_set
            pc_str = vals[13]  # e.g. struct('ACADS',2.44,'ACADVL',2.44)
            pc_str = "{" + pc_str.lstrip('struct(').rstrip(')')\
                .replace("',", "':") + "}"  # e.g. {'ACADS':2.44,'ACADVL':2.44}
            pway.gene_pc = eval(pc_str)
            allPathways.append(pway)
    return allPathways


def save_dissimilarity_files(pway_list_full, scores_path, names_path, min_len=50):
    """Reads gene_sets from all pathways, calculating a matrix of dissimilarity
    scores. Saves scores and path_id+names in text files at specified paths."""
    pway_list = [p for p in pway_list_full if p.gene_set]
    effect_sizes = [float(p.n_effective)/p.n_actual for p in pway_list]
    effect_sizes = misc.get_at_least_n(effect_sizes, min_len)  # top ~50 effects
    n_pways = len(effect_sizes)
    gene_set_list = [p.gene_set for p in pway_list][:n_pways]
    scores = misc.get_distance_matrix(gene_set_list)
    # SAVE SCORES SQUARE MATRIX IN TEXT FILE. (savetxt from numpy)
    np.savetxt(scores_path, scores, delimiter='\t')
    # SAVE NAMES IN TEXT FILE.
    names = [(p.path_id, p.nice_name) for p in pway_list if p.gene_set]
    names = names[:n_pways]  # trim down to length n_pways
    with open(names_path, 'w') as out:
        for n in names:
            out.write(n[0] + '\t' + n[1] + '\n')


def make_js_file(allPathways, out_path):
    # BUILD JS FILE
    with open(out_path, 'w') as out:  # 'pathways_pvalues_{}.js'
        out.write("root_url = '{}';\n".format(allPathways[0].root_url))
        out.write("pwayList = [\n")
        for ind, pway in enumerate(allPathways):
            if pway.gene_set:
                out.write(pway.as_string_js())
                if ind < len(allPathways) - 1:
                    out.write(",")
                out.write("\n")
        out.write("];\n")


def make_readable_file(allPathways, out_path):
    """Create file suitable for archiving for user."""
    header_str = "pathway_name, pathway_id, url, p_value, n_effective, n_actual, " \
                 "genes_mutated_pc, gene_len_min_kbp, shortest_gene(s), " \
                 "gene_len_max_kbp, longest_gene(s), " \
                 "gene_len_avg_kbp, gene_len_variance"
    header_line = '\t'.join(header_str.split(', ')) + '\n'
    with open(out_path, 'w') as out:
        out.write(header_line)
        for pway in allPathways:
            if pway.gene_set:
                line_vals = list()
                line_vals.append(pway.nice_name)
                line_vals.append(pway.path_id)
                line_vals.append(pway.full_url)
                line_vals.append(pway.p_value)
                line_vals.append(pway.n_effective)
                line_vals.append(pway.n_actual)
                line_vals.append(pway.gene_pc)
                line_vals.extend(list(pway.lengths_tuple))
                out.write('\t'.join([str(v) for v in line_vals]) + '\n')


def create_pway_plots(txt_path, scores_path, names_path, tree_path, genome_size,
                      hyper_path):
    """Run matlab script that builds matrix and target svgs."""
    # ORIG cmd = """matlab -nosplash -nodesktop -r "plot_pway_targets('{txtpath}');" < /dev/null >{root_dir}tempstdout.txt 2>{root_dir}tempstderr.txt &"""

    # worker = matlab.engine.start_matlab()
    plot_fn = "pway_plots({txtpath!r}, {scores!r}, {names!r}, {tree!r}, {gsize}, {hyper_path!r});".\
        format(txtpath=txt_path, scores=scores_path, names=names_path,
               tree=tree_path, gsize=int(genome_size), hyper_path=hyper_path)

    cmd = [current_app.config['MATLAB_PATH'], '-nosplash', '-nodesktop', '-r',
           plot_fn + ' exit;']

    with open(os.devnull, "r") as fnullin:
        with open(os.devnull, "w") as fnullout:
            proc = subprocess.Popen(cmd, stdin=fnullin,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
            (output, error) = proc.communicate()
            if error:
                raise MatlabFailureException('pway plot error: ' + error)
