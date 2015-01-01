"""Runs pathway pipeline on CancerDB tables or TCGA tables."""


from .decorators import async
from flask import current_app
from flask_login import current_user

import MySQLdb as mdb
from scipy.misc import comb
from numpy import *
from scipy import stats
import timeit
import argparse
import warnings
from collections import OrderedDict
import os
import subprocess

from . import db


dbvars = dict()

@async
def run_analysis_async(app, proj_dir, data_path, user_upload):
    """Asynchronous run of pathway analysis."""
    global dbvars
    with app.app_context():
        dbvars = dict(host=current_app.config['SGG_DB_HOST'],
                      db=current_app.config['SGG_DB_NAME'],
                      read_default_file=current_app.config['SGG_DB_CNF'])
        table_name = 'mutations_{}'.format(user_upload.file_id)
        MutationTable(table_name, data_path)
        run(proj_dir, table_name, user_upload)
        user_upload.run_complete = True
        db.session.add(user_upload)
        db.session.commit()


def run_analysis(proj_dir, data_path, user_upload):
    user_upload.run_accepted = True
    db.session.add(user_upload)
    db.session.commit()
    app = current_app._get_current_object()
    run_analysis_async(app, proj_dir, data_path, user_upload)


class NonSingleResult(Exception):
    pass


class Patient():
    def __init__(self, patient_id, n_mutated, is_mutated):
        self.patient_id = patient_id
        self.n_mutated = n_mutated
        self.is_mutated = is_mutated


class MutationTable():
    """Loads data from file into table, given table_name and file path."""
    def __init__(self, table_name, data_path):
        self.table_name = table_name
        self.data_path = data_path
        self.loaded = False
        # CREATE AND POPULATE TABLE
        create_str = """CREATE TABLE `{}` (
          `hugo_symbol` VARCHAR(255) DEFAULT NULL,
          `entrez_id` INT(11) DEFAULT NULL,
          `patient_id` VARCHAR(255) DEFAULT NULL,
          `variant_classification` VARCHAR(255) DEFAULT NULL,
          KEY `temp_patient_id` (`patient_id`) USING HASH,
          KEY `temp_patient_entrez` (`patient_id`,`entrez_id`)
          USING HASH);""".format(table_name)
        load_str = """load data local infile '{}'
        into table `{}` fields terminated by '\t'
        lines terminated by '\n' ignore 1 lines;""".format(
            data_path, table_name
        )
        try:
            con = mdb.connect(**dbvars)
            cur = con.cursor()
            cur.execute(create_str)
            cur.execute(load_str)
            con.commit()
            self.loaded = True
        except mdb.Error as e:
            print "Error %d: %s" % (e.args[0], e.args[1])
        finally:
            if con:
                con.close()



class GeneMatrix():
    """Holds patient-gene matrix info for a single pathway.
    Write matrix to file with call to export_matrix."""

    def __init__(self):
        self.genePatientDict = dict()

    def add_gene_patients(self, gene, patientList):
        """Add gene-patientList pairs into dictionary."""
        if gene in self.genePatientDict:
            self.genePatientDict.extend(patientList)
        else:
            self.genePatientDict[gene] = patientList

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
                 patient_ids=list(), max_mutations=500, expressed_table=None,
                 ignore_genes=list()):
        PathwaySummaryBasic.__init__(self, pathway_number)
        # above gives path_id, n_actual, n_effective, p_value
        self.patients = list()  # tuples. (patient_id, n_mutations, is_mutated)
        self.filter_patient_ids = patient_ids
        self.filter_expressed = expressed_table
        self.ignore_genes = ignore_genes
        self.max_mutations = max_mutations
        self.proj_abbrvs = proj_abbrvs
        self.gene_coverage = OrderedDict()
        self.exclusive_genes = list()
        self.cooccurring_genes = list()
        self.geneMatrix = None
        self.runtime = None  # only set up by file reader

    def set_up_from_file(self, pval, psize, peffect, runtime):
        """Manually specify pathway summary properties (after running init)."""
        self.p_value = pval
        self.n_actual = psize
        self.n_effective = peffect
        self.runtime = runtime
        # # fetch gene_coverage, exclusive_genes, cooccurring genes
        # self._populate_exclusive_cooccurring()
        # self._update_gene_coverage()

    def set_pathway_size(self):
        """Lookup pathway size in DB. Save to data attribute."""
        if self.ignore_genes:
            gene_filter = self._build_ignore_gene_filter_str(
                form="AND").replace("hugo_symbol", "symbol")
        else:
            gene_filter = ""
        pway_size = None
        # GET pway_size
        cmd1 = """SELECT count(DISTINCT pgl.entrez_id) AS pway_size 
        FROM refs.pathway_gene_link pgl
        INNER JOIN refs.ncbi_entrez n ON pgl.entrez_id = n.geneId
        {expression_filter_pgl}
        WHERE pgl.entrez_id IS NOT NULL {ignore_gene_filter}
        AND path_id = {pathid}
        GROUP BY path_id;""".format(pathid=self.path_id,
                                    expression_filter_pgl=
                                    self._build_expressed_filter_str(
                                        'pgl.entrez_id'),
                                    ignore_gene_filter=gene_filter)
        try:
            con = mdb.connect(**dbvars)
            cur = con.cursor()
            cur.execute(cmd1)
            # assert isinstance(cur.rowcount, int)
            row_count = cur.rowcount
            if not row_count == 1:
                print "Result contains {} rows Ids for pathway {}.".format(
                    row_count, self.path_id)
                self.n_actual = 0
                return
            pway_size = cur.fetchone()
        except mdb.Error as e:
            print "Error %d: %s" % (e.args[0], e.args[1])
        finally:
            if con:
                con.close()
        self.n_actual = int(pway_size[0])

    def populate_patient_info(self):
        """Fetch info on all patients from DB, add to patients list."""
        if self.proj_abbrvs:
            self._append_patients()

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

    def _append_patients(self):
        """Gets patient-pathway gene overlap info from database.
        # Returns tuple collection of row tuples.
        Row tuple: (PATIENT_ID, N_PATIENT, BOOL_MUTATED), e.g. (678L, 323L, 1L)
        """
        # in case of genes in ignore list, choose prefix for sql snippet
        if self.filter_patient_ids:
            ignore_gene_prefix = "AND"
        else:
            ignore_gene_prefix = "WHERE"
        if self.ignore_genes or self.filter_patient_ids:
            variant_where_and = "AND"
        else:
            variant_where_and = "WHERE"
        for tcga_table in self.proj_abbrvs:
            rows = None
            # # GET mutated boolean for genes in specified pathway for GOOD
            # patients (<max_mutations), even those without somatic mutation.
            cmd2 = """SELECT p.patient_id, n_patient,
                g.patient_id IS NOT NULL AS mutated FROM
                # good patients, mutation_count
                (SELECT patient_id, count(DISTINCT m.entrez_id) AS n_patient
                FROM {table} m
                    {expression_filter_m}
                    {patient_filter} {ignore_gene_filter1}
                    {variant_WHERE_AND} `Variant_Classification` <> 'Silent'
                    GROUP BY patient_id HAVING count(*) <= {max_mutations}) p
                LEFT JOIN
                # pathway mutation counts above 0
                (SELECT patient_id FROM
                    # get distinct genes from patients
                    (SELECT DISTINCT patient_id, entrez_id FROM {table}
                        {patient_filter} {ignore_gene_filter2}) pg
                    # join to pathway genes
                        INNER JOIN 
                    (SELECT pgl.entrez_id FROM refs.pathway_gene_link pgl 
                        {expression_filter_pgl}
                        WHERE path_id = {path_id}) pwg 
                    ON pwg.entrez_id = pg.entrez_id 
                    GROUP BY `patient_id`) g
                ON p.patient_id = g.patient_id;""" \
                .format(path_id=self.path_id,
                        max_mutations=self.max_mutations,
                        table=tcga_table,
                        patient_filter=self._build_patient_filter_str(
                            form='WHERE'),
                        ignore_gene_filter1=self._build_ignore_gene_filter_str(
                            form=ignore_gene_prefix),
                        ignore_gene_filter2=self._build_ignore_gene_filter_str(
                            form=ignore_gene_prefix),
                        expression_filter_m=self._build_expressed_filter_str(
                            'm.entrez_id'),
                        expression_filter_pgl=self._build_expressed_filter_str(
                            'pgl.entrez_id'),
                        variant_WHERE_AND=variant_where_and)
            try:
                con = mdb.connect(**dbvars)
                cur = con.cursor()
                cur.execute(cmd2)
                rows = cur.fetchall()

            except mdb.Error as e:
                print "Error %d: %s" % (e.args[0], e.args[1])
            finally:
                if con:
                    con.close()

            for row in rows:
                patient_id = row[0]
                n_mutated = row[1]
                is_mutated = row[2] > 0
                is_yale = False
                patient = Patient(patient_id, n_mutated, is_mutated)
                self.patients.append(patient)

    def update_exclusive_cooccurring_coverage(self):
        """POSTPROCESSING.
        Gather pathway gene info and write detailed output."""
        self.geneMatrix = GeneMatrix()  # populated during update_gene_coverage
        self._populate_exclusive_cooccurring()
        self._update_gene_coverage()

    def _populate_exclusive_cooccurring(self):
        """ Postprocessing step. Look up gene combinations hit
        (via _get_gene_combs_hit_yale and _get_gene_combs_hit_tcga),
        sort into sets."""
        gene_combs_list = self._get_gene_combs_hit()
        all_hit_genes_set = set()
        exclusive_gene_set = set()
        # get minimum set of gene combinations
        temp_set = set()
        for genes in gene_combs_list:
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

    def _get_gene_combs_hit(self):
        """Gets patient-pathway gene overlap info from databse.
        Only called by _populate_exclusive_cooccurring.
        """
        rows = None
        geneLists = list()
        for tcga_table in self.proj_abbrvs:
            # UNIQUE HUGO LISTS
            cmd = """SELECT DISTINCT symbols FROM
                (
                # PATIENT, HUGO PAIRS in pathway of interest.
                SELECT patient_id, group_concat(DISTINCT hugo_symbol
                ORDER BY hugo_symbol SEPARATOR ',') AS symbols
                FROM {table} t
                # gene subset in pathway of interest
                INNER JOIN refs.`pathway_gene_link` pgl
                ON t.entrez_id = pgl.entrez_id {expression_filter_pgl}
                NATURAL JOIN
                #good patients
                (SELECT patient_id FROM {table}
                WHERE `Variant_Classification` <> 'Silent' {patient_filter}
                    GROUP BY patient_id 
                    HAVING count(*) <= {max_mutations}) p
                WHERE path_id = {path_id} AND `Variant_Classification`<>'Silent'
                GROUP BY patient_id
                ) g;""". \
                format(path_id=self.path_id,
                       max_mutations=self.max_mutations,
                       table=tcga_table,
                       patient_filter=self._build_patient_filter_str(
                           form="AND"),
                       expression_filter_pgl=self._build_expressed_filter_str(
                           "pgl.entrez_id"))
            try:
                con = mdb.connect(**dbvars)
                cur = con.cursor()
                cur.execute(cmd)
                rows = cur.fetchall()
            except mdb.Error as e:
                print "Error %d: %s" % (e.args[0], e.args[1])
            finally:
                if con:
                    con.close()
            gene_strings = [row[0] for row in rows]
            temp_gene_lists = [genesString.split(',') for genesString in
                               gene_strings]
            geneLists.extend(temp_gene_lists)
        return geneLists

    def _update_gene_coverage(self):
        """Populate object's gene_coverage dictionary."""
        gene_coverage = dict()
        (total_patients, counts_temp) = self._get_gene_counts()
        for gene in counts_temp.iterkeys():
            if gene in gene_coverage:
                gene_coverage[gene] += counts_temp[gene]
            else:
                gene_coverage[gene] = counts_temp[gene]

        coverage_tuple = [tuple([gene, gene_coverage[gene]]) for gene in
                          gene_coverage]
        coverage_tuple.sort(key=lambda x: x[0])  # alphabeticize
        coverage_tuple.sort(key=lambda x: x[1],
                            reverse=True)  # sort by n_patients
        gene_coverage = OrderedDict(coverage_tuple)

        # convert to coverage by dividing by number of patients
        for gene in gene_coverage:
            gene_coverage[gene] = float(
                gene_coverage[gene]) / total_patients * 100
        self.gene_coverage = gene_coverage
        return

    def _get_gene_counts(self):
        """ Fetch dictionary: gene -> int percentage of patients with mutation.
        Dictionary may be empty if no pathway genes were mutated."""
        total_patients = 0  # will increase as inspect tables
        count_dict = dict()
        rows = None
        for tcga_table in self.proj_abbrvs:
            cmd0 = """SET SESSION group_concat_max_len = 30000;"""
            cmd1 = """SELECT count(DISTINCT patient_id) AS n_patient FROM 
                (SELECT patient_id FROM {table}
                WHERE `Variant_Classification` <> 'Silent' {patient_filter}
                GROUP BY patient_id HAVING count(*) <= {max_mutations}) p2;""" \
                .format(max_mutations=self.max_mutations, table=tcga_table,
                        patient_filter=self._build_patient_filter_str(
                            form="AND"))
            # HUGO LIST AND PATIENT COUNTS
            cmd2 = """SELECT hugo_symbol, count(DISTINCT patient_id)
                AS n_patients, GROUP_CONCAT(DISTINCT patient_id) AS patients
                FROM {table} t
                # gene subset in pathway of interest
                INNER JOIN refs.`pathway_gene_link` pgl
                ON t.entrez_id = pgl.entrez_id {expression_filter_pgl}
                NATURAL JOIN
                #good patients
                (SELECT patient_id FROM {table} t
                    WHERE `Variant_Classification` <> 'Silent' {patient_filter}
                    GROUP BY patient_id HAVING count(*) <= {max_mutations} )  p
                WHERE path_id = {path_id}
                    AND `Variant_Classification` <> 'Silent'
                {patient_filter}
                GROUP BY hugo_symbol
                ORDER BY hugo_symbol;""" \
                .format(path_id=self.path_id,
                        max_mutations=self.max_mutations,
                        table=tcga_table,
                        patient_filter=self._build_patient_filter_str(
                            form="AND"),
                        expression_filter_pgl=self._build_expressed_filter_str(
                            "pgl.entrez_id"))
            try:
                con = mdb.connect(**dbvars)
                cur = con.cursor()
                cur.execute(cmd0)
                cur.execute(cmd1)
                rowCount = cur.rowcount
                if not rowCount == 1:
                    raise NonSingleResult(
                        "Result contains %g rows Ids for pathway %s."
                        % (rowCount, self.path_id))
                total_patients += int(cur.fetchone()[0])
                cur.execute(cmd2)
                rowCount = cur.rowcount
                if not rowCount:
                    # NO GENES MUTATED. n_effective < n_pathway
                    return total_patients, count_dict
                rows = cur.fetchall()
            except mdb.Error as e:
                print "Error %d: %s" % (e.args[0], e.args[1])
            finally:
                if con:
                    con.close()
            for geneRow in rows:
                gene = geneRow[0]
                coverage = int(geneRow[1])
                patient_names = geneRow[2].split(',')
                if not len(patient_names) == coverage:
                    raise Exception(
                        "Pathway coverage query gives inconsistent " +
                        "patient counts and patient names.")
                self.geneMatrix.add_gene_patients(gene, patient_names)
                if gene in count_dict:
                    count_dict[gene] += coverage
                else:
                    count_dict[gene] = coverage
        return total_patients, count_dict


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
        self.max_mutations = pway.max_mutations

    def run(self):
        """ Calculate likelihood and maximum likelihood estimate."""
        self.likelihood = self._get_pway_likelihood(self.pway.n_actual)
        (ne, lastll) = self._get_ne()
        self.ne = ne
        self.ne_ll = lastll
        self.D = -2 * self.likelihood + 2 * self.ne_ll
        self.pvalue = 1 - stats.chi2.cdf(self.D, 1)
        self.pway.n_effective = self.ne
        self.pway.p_value = self.pvalue

    def _get_pway_likelihood(self, pway_size=None):
        if pway_size is None:
            pway_size = self.pway_size
        prob_list = list()
        # iterate over patients
        for patient in self.pway.patients:
            n_patient = patient.n_mutated
            mutated = patient.is_mutated
            # if there are fewer out-of-pathway genes than there are
            # genes mutated, p_no_mut is zero
            if self.G - pway_size < n_patient:
                p_no_mut = float64(0)
            else:
                p_no_mut = exp(math.log(comb(self.G - pway_size, n_patient,
                                             exact=True)) - math.log(
                    comb(self.G, n_patient, exact=True)))
            if mutated:
                p = 1 - p_no_mut
            else:
                p = p_no_mut
            prob_list.append(log(p))
        prob_array = array(prob_list)
        return prob_array.sum()

    def _get_ne(self):
        last_ll = None
        # improved = False
        ne = None
        # profile = list()
        # if pathway_size is zero, effective size is zero.
        if not self.pway.n_actual:
            ne = 0
            ult = float64(0)
            warnings.warn(
                "Pathway {} contains zero genes. ".format(self.pway.path_id))
            return ne, ult
        # if all patients mutated, use ne=genome_size - max_mutations
        if False not in [patient.is_mutated for patient in self.pway.patients]:
            ne = self.G
            ult = float64(0)
            warnings.warn("All patients have mutation in pathway {}. ".format(
                self.pway.path_id) + "Effective size is full genome.")
            return (ne, ult)
        # check last 2 vals to check for decline:
        penult = self._get_pway_likelihood(
            pway_size=self.G - 2)  # WAS self.G - self.max_mutations - 1
        ult = self._get_pway_likelihood(
            pway_size=self.G - 1)  # WAS self.G - self.max_mutations
        if ult > penult:
            ne = self.G
            return (ne, ult)
        # at this stage, there will be a max before Genome size
        for pway_size in xrange(1,self.G):
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
                    warnings.warn("Premature stop for pway {}.".format(
                        self.pway.path_id))
                break
            last_ll = this_ll
        return (ne, last_ll)


class GenericPathwayFileProcessor():
    """Generic object that can convert yale_proj_ids and tcga_proj_abbrvs
    to file_name."""
    base_str = 'pathways_pvalues_'

    def __init__(self, dir_path, file_id, name_suffix=None):
        self.file_id = file_id
        self.dir_path = dir_path
        self.name_suffix = name_suffix
        self.root_name = self._get_root_filename()

    def _get_root_filename(self):
        """Root file name for output files."""
        root_name = os.path.join(self.dir_path, self.base_str
                                 + str(self.file_id))
        if self.name_suffix:
            root_name += '_' + self.name_suffix
        return root_name


class PathwayBasicFileWriter(GenericPathwayFileProcessor):
    """Writes initial p-value file."""

    def write_pvalue_file(self, lcalc, runtime):
        """Write initial processing file with p-value and MLE estimate."""
        bufsize = 1  # line buffered output
        outfile_name = self.root_name + '.txt'
        path_id = lcalc.pway.path_id
        with open(outfile_name, 'a', bufsize) as out:
            out.write('{}\t{:.3e}\t{}\t{}\t{:.2f}\n'.format(
                path_id, lcalc.pvalue, lcalc.pway.n_actual, lcalc.ne, runtime))


class PathwayListAssembler(GenericPathwayFileProcessor):
    """Builds ordered list of pathways from basic p-value file."""

    def __init__(self, dir_path, file_id, proj_abbrvs, patient_ids=list(),
                 name_suffix=None, max_mutations=None, expressed_table=None,
                 ignore_genes=list()):
        # create self.root_name
        GenericPathwayFileProcessor.__init__(self, dir_path, file_id,
                                             name_suffix=name_suffix)
        self.proj_abbrvs = proj_abbrvs
        self.filter_patient_ids = patient_ids
        self.max_mutations = max_mutations
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
                runtime = float(row[4])
                # set up pathway object
                pway = PathwaySummary(path_id, self.proj_abbrvs,
                                      patient_ids=self.filter_patient_ids,
                                      max_mutations=self.max_mutations,
                                      expressed_table=self.expressed_table,
                                      ignore_genes=self.ignore_genes)
                pway.set_up_from_file(pval, psize, peffect, runtime)
                all_pathways.append(pway)
        # sort pathways by effect_size : actual_size
        all_pathways.sort(
            key=lambda pway: self._get_size_ratio(pway.n_effective,
                                                  pway.n_actual),
            reverse=True)
        # sort pathways by p-value
        all_pathways.sort(key=lambda pw: pw.p_value, reverse=False)
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
                 name_suffix=None):
        # create self.root_name
        GenericPathwayFileProcessor.__init__(self, dir_path, file_id,
                                             name_suffix=name_suffix)
        self.allPathways = pway_object_list
        self.nameDict = self.get_pathway_name_dict()
        self.outfile_name = self.root_name + self.name_postfix
        self.matrix_folder = 'matrix_txt'

    @staticmethod
    def get_pathway_name_dict():
        """Gets name for all pathways, stored in dict: pathid -> pathname."""
        rows = None
        pathway_dict = dict()
        # GET pway_size
        cmd1 = """SELECT p.path_id, pathway_name FROM refs.pathways p
        INNER JOIN
        (SELECT DISTINCT path_id FROM refs.pathway_gene_link) l
        ON p.path_id = l.path_id;"""
        try:
            con = mdb.connect(**dbvars)
            cur = con.cursor()
            cur.execute(cmd1)
            row_count = cur.rowcount
            if not row_count > 1:
                raise Exception(
                    "Result contains %g rows Ids for pathway lookup."
                    % row_count)
            rows = cur.fetchall()
        except mdb.Error as e:
            print "Error %d: %s" % (e.args[0], e.args[1])
        finally:
            if con:
                con.close()
        # rows is [[id,name],[id,name],...]
        for pair in rows:
            path_id = int(pair[0])
            path_name = pair[1]
            pathway_dict[path_id] = path_name
        return pathway_dict

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
                if pway.p_value < 0.1:
                    pway.update_exclusive_cooccurring_coverage()
                    if pway.gene_coverage:  # if genes are hit...
                        if pway.n_effective > pway.n_actual:
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
                out.write("""{path_id}\t{name}\t{n_actual}\t{n_effective}\t{p_value:.3e}\t{runtime:.2f}\t{exclusive_string}\t{cooccurring_string}\t{coverage_string}\n"""
                          .format(path_id=pway.path_id, name=path_name,
                                  coverage_string=coverage_string,
                                  n_actual=pway.n_actual,
                                  n_effective=pway.n_effective,
                                  p_value=pway.p_value,
                                  runtime=pway.runtime,
                                  exclusive_string=exclusive_string,
                                  cooccurring_string=cooccurring_string))

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
        self.name = None
        self.root_url = "http://www.broadinstitute.org/gsea/msigdb/cards/"
        self.url = self.fetch_url()
        self.gene_set = set()

    def as_string(self):
        """Return string for javascript."""
        name_str = self.name.replace('_', ' ')
        name_str = name_str.replace('BIOCARTA ', '')
        name_str = name_str.replace('REACTOME ', '')
        name_str = name_str.replace('KEGG ', '')
        name_str = name_str.replace('TEL PATHWAY', 'TELOMERASE PATHWAY')
        name_str = name_str.replace('RNA PATHWAY', 'PKR SIGNALING PATHWAY')
        outstr = "{ id:" + self.path_id + ", " + "name:'" + name_str \
                 + "', pval:'" + self.p_value + "', size:" + str(self.n_actual) \
                 + ", effective:" + str(
            self.n_effective) + ", url:'" + self.url + "', geneSet: "
        if self.gene_set:
            gene_set_str = "','".join(self.gene_set)
            outstr = outstr + "['" + gene_set_str + "']}"
        else:
            outstr = outstr + "[]}"
        return outstr

    def fetch_url(self):
        """Look up database (refs.pathways) to get info_url for path_id."""
        url_row = None
        cmd = """SELECT info_url FROM refs.pathways WHERE path_id = {};""".format(
            self.path_id)
        try:
            con = mdb.connect(**dbvars)
            cur = con.cursor()
            cur.execute(cmd)
            rowCount = cur.rowcount
            if not rowCount == 1:
                raise Exception("Result contains %g rows Ids for pathway %s."
                                % (rowCount, self.path_id))
            url_row = cur.fetchone()
        except mdb.Error as e:
            print "Error %d: %s" % (e.args[0], e.args[1])
        finally:
            if con:
                con.close()
        url = url_row[0]
        url = url.split(self.root_url)[1]
        return url


class PathwayIdsFetcher():
    def __init__(self, interest_genes):
        self.interest_genes = interest_genes

    def fetchPathwayIds(self):
        """Get pathway ids containing genes in (possibly empty) interest set."""
        rows = None
        all_path_ids = list()
        if self.interest_genes:
            genes_string = repr(tuple(self.interest_genes))
            genes_string = genes_string.replace(",)", ")")
            genes_string = "WHERE symbol IN " + genes_string
        else:
            genes_string = ""
        # GET pway_size
        cmd1 = """SELECT distinct path_id FROM refs.pathway_gene_link pgl
            INNER JOIN (SELECT geneId FROM refs.ncbi_entrez 
                {genes_string}) g
            ON pgl.entrez_id = g.geneId ORDER BY path_id;""".format(
            genes_string=genes_string)
        try:
            con = mdb.connect(**dbvars)
            cur = con.cursor()
            cur.execute(cmd1)
            rowCount = cur.rowcount
            if not rowCount:
                raise Exception(
                    "Result contains %g rows Ids for pathway lookup."
                    % (rowCount))
            rows = cur.fetchall()
        except mdb.Error as e:
            print "Error %d: %s" % (e.args[0], e.args[1])
        finally:
            if con:
                con.close()
        # rows is [[id,name],[id,name],...]
        for id in rows:
            all_path_ids.append(int(id[0]))
        return all_path_ids


class BackgroundGenomeFetcher():
    def __init__(self, genome_str, expressed_table=None):
        """Specifying an expressed table will result in a genome size equal to 
        the number of expressed genes, unless a genome_str is specified. The 
        default genome_str for runs without expression data is no_pseudo."""
        if not genome_str and not expressed_table:
            genome_str = 'no_pseudo'
        self.genome_size = self._fetch_genome_size(genome_str, expressed_table)

    def _fetch_genome_size(self, genome_str, expressed_table):
        # genome_size
        if genome_str == 'protein-coding':
            genome_size = 20462
        elif genome_str == 'no_pseudo':
            genome_size = 28795
        elif genome_str == 'all':
            genome_size = 45466
        elif genome_str == 'inc_misc_chr':
            genome_size = 46286
        elif expressed_table:
            genome_size = self._fetch_expressed_genome_size(expressed_table)
        else:
            raise Exception("Unknown genome version")
        print(
            "Using genome version '{}': {} genes".format(genome_str,
                                                         genome_size))
        return genome_size

    @staticmethod
    def _fetch_expressed_genome_size(expressed_table):
        """Count genes via SQL query: assumes row count equals gene count."""
        cmd1 = """SELECT count(*) FROM {table_name};""".format(
            table_name=expressed_table)
        try:
            con = mdb.connect(**dbvars)
            cur = con.cursor()
            cur.execute(cmd1)
            rowCount = cur.rowcount
            if not rowCount or rowCount > 1:
                raise Exception("Expressed genome size db-lookup failed.")
            rows = cur.fetchall()
        except mdb.Error as e:
            print "Error %d: %s" % (e.args[0], e.args[1])
        finally:
            if con:
                con.close()
        return int(rows[0][0])







def run(dir_path, table_name, user_upload):
    """ EXAMPLE ARGUMENTS
    dir_path = /Users/sgg/Downloads/uploads/1/41/
    table_name = mutations_41
    user_upload is upload object (file_id, user_id, filename, ...etc)
    """

    # max_mutations
    file_id = user_upload.file_id
    table_list = [table_name]  # code can iterate through list of tables
    max_mutations = user_upload.n_cutoff or 500  # default max is 500
    if user_upload.ignore_genes:
        ignore_genes = user_upload.ignore_genes.split(',')
    else:
        ignore_genes = []
    genome_size = BackgroundGenomeFetcher(user_upload.genome_size,
                                          None).genome_size
    proj_suffix = user_upload.proj_suffix
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
        interest_genes = tuple(user_upload.required_genes.split(','))
    else:
        interest_genes = tuple()
    print("Interest genes: {}".format(str(interest_genes)))

    id_fetcher = PathwayIdsFetcher(interest_genes)
    all_path_ids = id_fetcher.fetchPathwayIds()

    for pathway_number in all_path_ids[0:30]:
        # Populate pathway object, and time pvalue calculation
        start = timeit.default_timer()
        pway = PathwaySummary(pathway_number, table_list,
                              max_mutations=max_mutations,
                              expressed_table=None,
                              ignore_genes=ignore_genes)
        pway.set_pathway_size()
        pway.populate_patient_info()
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
                                     max_mutations=max_mutations,
                                     # expressed_table=args.expression,
                                     ignore_genes=ignore_genes)
    pway_list = assembler.get_ordered_pway_list()

    # Rank pathways, gather extra stats and write to final file
    final_writer = PathwayDetailedFileWriter(dir_path, file_id, pway_list,
                                             name_suffix=proj_suffix)
    final_writer.write_detailed_file()

    # generate_files(dir_path, final_writer.outfile_name, user_upload)

    detail_path = final_writer.outfile_name
    descriptive_name = str(user_upload.file_id) + user_upload.filename
    js_out_path = os.path.join(dir_path, descriptive_name + ".js")

    make_js_file(detail_path, js_out_path)
    # html_name = create_html(detail_path, out_dir, project_str, descriptive_name,
    #                         skipfew)
    create_svgs(detail_path)
    create_matrix_svgs(detail_path)




def make_js_file(results_path, out_path):
    # BUILD JS FILE
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
            if (skipLine):
                continue
            pway.name = name
            pway.n_actual = int(vals[2])
            pway.n_effective = int(vals[3])
            pway.p_value = vals[4]
            gene_set = set.union(set(eval(vals[6])), set(eval(vals[7])))
            pway.gene_set = gene_set
            allPathways.append(pway)
    with open(out_path, 'w') as out:  # 'pathways_pvalues_{}.js'
        out.write("root_url = '{}';\n".format(pway.root_url))
        out.write("pwayList = [\n")
        for ind, pway in enumerate(allPathways):
            if pway.gene_set:
                out.write(pway.as_string())
                if ind < len(allPathways) - 1:
                    out.write(",")
                out.write("\n")
        out.write("];\n")


def create_svgs(txt_path):
    """Run matlab script that builds svgs in directory containing
    pathways txt."""
    # ORIG cmd = """matlab -nosplash -nodesktop -r "plot_pway_targets('{txtpath}');" < /dev/null >{root_dir}tempstdout.txt 2>{root_dir}tempstderr.txt &"""
    cmd = 'matlab -nosplash -nodesktop -r \"plot_pway_targets({txtpath!r},' \
          '\'--svg\',\'--skipfew\');\"'.format(txtpath=txt_path)
    with open(os.devnull, "r") as fnullin:
        with open(os.devnull, "w") as fnullout:
            subprocess.check_call(cmd, stdin=fnullin, stdout=fnullout,
                                  stderr=fnullout, shell=True)


def create_matrix_svgs(txt_path):
    """Run matlab script that builds matrix svgs. SVGs stored in matrix_
    txt/matrix_svg."""
    # ORIG cmd = """matlab -nosplash -nodesktop -r "plot_pway_targets('{txtpath}');" < /dev/null >{root_dir}tempstdout.txt 2>{root_dir}tempstderr.txt &"""
    cmd = 'matlab -nosplash -nodesktop -r \"plot_patient_genes(' \
          '{txtpath!r});\"'.format(txtpath=txt_path)
    with open(os.devnull, "r") as fnullin:
        with open(os.devnull, "w") as fnullout:
            subprocess.check_call(cmd, stdin=fnullin, stdout=fnullout,
                                  stderr=fnullout, shell=True)


if __name__ == '__main__':
    main()

    


