#!/Users/Stephen/Library/Enthought/Canopy_64bit/User/bin/python

"""Runs pathway pipeline on CancerDB tables or TCGA tables."""

import MySQLdb as mdb
from scipy.misc import comb
from numpy import *
from scipy import stats
import timeit
import argparse
import warnings
from collections import OrderedDict
import os


class NonSingleResult(Exception):
    pass


class Patient():
    def __init__(self, patient_id, n_mutated, is_mutated, is_yale):
        self.patient_id = patient_id
        self.n_mutated = n_mutated
        self.is_mutated = is_mutated
        self.is_yale = is_yale


class GeneMatrix():
    """Holds patient-gene matrix info for a single pathway."""

    def __init__(self):
        self.genePatientDict = dict()

    def addGenePatients(self, gene, patientList):
        """Add gene-patientList pairs into dictionary."""
        if gene in self.genePatientDict:
            self.genePatientDict.extend(patientList)
        else:
            self.genePatientDict[gene] = patientList

    def export_matrix(self, outfile, exclusive_genes):
        """Writes tab-separated matrix file for patient/gene pairs in pathway."""
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


class PathwaySummary():
    """Holds pathway information, and can fetch info from db."""

    def __init__(self, pathway_number, yale_proj_ids, tcga_proj_abbrvs,
                 patient_ids=list(), max_mutations=500, expressed_table=None,
                 ignore_genes=list()):
        self.path_id = pathway_number
        self.n_actual = None
        self.patients = list()  # tuples. (patient_id, n_mutations, is_mutated)
        self.filter_patient_ids = patient_ids
        self.filter_expressed = expressed_table
        self.ignore_genes = ignore_genes
        self.n_effective = None
        self.p_value = None
        self.max_mutations = max_mutations
        self.yale_proj_ids = yale_proj_ids
        self.tcga_proj_abbrvs = tcga_proj_abbrvs
        # populated during post-processing:
        self.gene_coverage = OrderedDict()
        self.exclusive_genes = list()
        self.cooccurring_genes = list()
        self.geneMatrix = None  # set during _populate_exclusive_cooccurring_coverage()
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
        if self.yale_proj_ids:
            self._append_yale_patients()
        if self.tcga_proj_abbrvs:
            self._append_tcga_patients()

    def _build_patient_filter_str(self, form="WHERE"):
        """build SQL substring to filter patients by ids in mutation lookup."""
        if form == 'WHERE':
            prepend = "WHERE "
        elif form == 'AND':
            prepend = "AND "
        else:
            raise Exception("Unrecognized form in patient filter.")
        if self.filter_patient_ids:
            filter = (prepend + "patient_id IN " +
                      str(self.filter_patient_ids).replace("[", "(").replace(
                          "]", ")"))
        else:
            filter = ""
        return filter

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
        join_ref is abbreviation of table and column to join to expression table,
        e.g. m.entrez_id or pwg.entrez_gene_id.
        Assumes expression table is in tcga database.
        """
        if self.filter_expressed:
            filter = "INNER JOIN tcga.{filter_expressed} e ON {join_ref} = e.entrez_id" \
                .format(filter_expressed=self.filter_expressed,
                        join_ref=join_ref)
        else:
            filter = ""
        return filter

    def _append_yale_patients(self):
        """Get patient-pathway gene overlap info from database.
        # Returns tuple collection of row tuples.
        Row tuple: (PATIENT_ID, N_PATIENT, BOOL_MUTATED), e.g. (678L, 323L, 1L)
        """
        rows = None
        # # GET mutated boolean for genes in specified pathway for GOOD
        # patients (<max_mutations), even those without somatic mutation.
        cmd2 = """SELECT p.patient_id, n_patient,
            g.patient_id IS NOT NULL AS mutated FROM
            # good patients for project, mutation_count
            (SELECT patient_id, count(DISTINCT m.entrez_gene_id) AS n_patient
            FROM
                mutations_tumor_normal m NATURAL JOIN normals
                NATURAL JOIN patients {expression_filter_m}
                WHERE project_id IN ({projGroupStr})
                AND `Variant_Classification` <> 'Silent'
                {patient_filter} {ignore_gene_filter}
            GROUP BY patient_id HAVING count(*) <= {max_mutations})  p 
            LEFT JOIN
            # patients in project with mutation in pathway
            (SELECT patient_id FROM 
                    (SELECT entrez_id FROM refs.`pathway_gene_link` pwg 
                        WHERE path_id={path_id}) pwg 
                INNER JOIN mutations_tumor_normal m
                    ON pwg.entrez_id=m.entrez_gene_id
                NATURAL JOIN normals NATURAL JOIN patients
                {expression_filter_m} 
                WHERE project_id IN ({projGroupStr}
                AND `Variant_Classification` <> 'Silent' {ignore_gene_filter})
                GROUP BY patient_id) g 
            ON p.patient_id = g.patient_id;""" \
            .format(path_id=self.path_id,
                    max_mutations=self.max_mutations,
                    projGroupStr=','.join(
                        str(i) for i in
                        self.yale_proj_ids),
                    patient_filter=self._build_patient_filter_str(
                        form="AND"),
                    ignore_gene_filter=self._build_ignore_gene_filter_str(
                        form="AND"),
                    expression_filter_m=self._build_expressed_filter_str(
                        'm.entrez_gene_id'))
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
            is_yale = True
            patient = Patient(patient_id, n_mutated, is_mutated, is_yale)
            self.patients.append(patient)

    def _append_tcga_patients(self):
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
            variant_WHERE_AND = "AND"
        else:
            variant_WHERE_AND = "WHERE"
        for tcga_table in self.tcga_proj_abbrvs:
            rows = None
            # # GET mutated boolean for genes in specified pathway for GOOD
            # patients (<max_mutations), even those without somatic mutation.
            cmd2 = """SELECT p.patient_id, n_patient,
                g.patient_id IS NOT NULL AS mutated FROM
                # good patients, mutation_count
                (SELECT patient_id, count(DISTINCT m.entrez_id) AS n_patient
                FROM tcga.{table} m
                    {expression_filter_m}
                    {patient_filter} {ignore_gene_filter1}
                    {variant_WHERE_AND} `Variant_Classification` <> 'Silent'
                    GROUP BY patient_id HAVING count(*) <= {max_mutations}) p
                LEFT JOIN
                # pathway mutation counts above 0
                (SELECT patient_id FROM
                    # get distinct genes from patients
                    (SELECT DISTINCT patient_id, entrez_id FROM tcga.{table}
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
                        variant_WHERE_AND=variant_WHERE_AND)
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
                patient = Patient(patient_id, n_mutated, is_mutated, is_yale)
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
        gene_combs_list = list()
        all_hit_genes_set = set()
        exclusive_gene_set = set()
        if self.yale_proj_ids:
            tempList = self._get_gene_combs_hit_yale()
            gene_combs_list.extend(tempList)
        if self.tcga_proj_abbrvs:
            tempList = self._get_gene_combs_hit_tcga()
            gene_combs_list.extend(tempList)
        # get minimum set of gene combinations
        tempSet = set()
        for genes in gene_combs_list:
            tempSet.add(tuple(genes))
        gene_combs_list = [list(genes) for genes in tempSet]
        for geneList in gene_combs_list:
            if len(geneList) == 1:
                exclusive_gene_set.add(geneList[0])
            for gene in geneList:
                all_hit_genes_set.add(gene)
        cooccurring_gene_set = all_hit_genes_set.difference(exclusive_gene_set)
        # onlyPairedWithExclusive = dict(zip(cooccurring_gene_set,
        # [True for i in cooccurring_gene_set]))
        # for geneList in gene_combs_list:
        # if not exclusive_gene_set.intersection(set(geneList)):
        # for gene in geneList:
        # onlyPairedWithExclusive[gene] = False
        # tagAlongGenes = [gene for gene in cooccurring_gene_set 
        # if onlyPairedWithExclusive[gene]]
        # coExclusiveGenes = [gene for gene in cooccurring_gene_set if not 
        # onlyPairedWithExclusive[gene]]
        # for gene in coExclusiveGenes:
        # exclusive_gene_set.add(gene)
        self.exclusive_genes = sorted(list(exclusive_gene_set))
        # self.cooccurring_genes = sorted(tagAlongGenes)
        self.cooccurring_genes = sorted(list(cooccurring_gene_set))
        return

    def _get_gene_combs_hit_yale(self):
        """Gets patient-pathway gene overlap info from databse.
        Only called by _populate_exclusive_cooccurring.
        """
        rows = None
        projIdsTuple = self.yale_proj_ids

        # UNIQUE HUGO LISTS
        cmd = """SELECT DISTINCT symbols FROM
            (# PATIENT, HUGO PAIRS in pathway of interest.
            SELECT patient_id, group_concat(DISTINCT hugo_symbol
            ORDER BY hugo_symbol SEPARATOR ',') AS symbols
            FROM 
                # gene subset in pathway of interest            
                (SELECT pgl.entrez_id FROM refs.`pathway_gene_link` pgl
                    {expression_filter_pgl}
                    WHERE path_id = {path_id}) pgl
            INNER JOIN mutations_tumor_normal m 
            ON m.entrez_gene_id = pgl.entrez_id
            NATURAL JOIN normals 
            NATURAL JOIN patients p
            WHERE project_id IN ({projGroupStr})
            {patient_filter} AND `Variant_Classification` <> 'Silent'
            GROUP BY patient_id
            ) g
            INNER JOIN
            #good patients
            (SELECT patient_id FROM mutations_tumor_normal NATURAL JOIN normals
                NATURAL JOIN patients WHERE project_id IN ({projGroupStr})
                {patient_filter} AND `Variant_Classification` <> 'Silent'
                GROUP BY patient_id HAVING count(*) <= {max_mutations} )  p2
            ON g.patient_id = p2.patient_id;""". \
            format(path_id=self.path_id,
                   max_mutations=self.max_mutations,
                   projGroupStr=','.join(
                       str(i) for i in
                       projIdsTuple),
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
        geneStrings = [row[0] for row in rows]
        geneLists = [genesString.split(',') for genesString in geneStrings]
        return geneLists

    def _get_gene_combs_hit_tcga(self):
        """Gets patient-pathway gene overlap info from databse.
        Only called by _populate_exclusive_cooccurring.
        """
        rows = None
        geneLists = list()
        for tcga_table in self.tcga_proj_abbrvs:
            # UNIQUE HUGO LISTS
            cmd = """SELECT DISTINCT symbols FROM
                (
                # PATIENT, HUGO PAIRS in pathway of interest.
                SELECT patient_id, group_concat(DISTINCT hugo_symbol
                ORDER BY hugo_symbol SEPARATOR ',') AS symbols
                FROM tcga.{table} t
                # gene subset in pathway of interest
                INNER JOIN refs.`pathway_gene_link` pgl
                ON t.entrez_id = pgl.entrez_id {expression_filter_pgl}
                NATURAL JOIN
                #good patients
                (SELECT patient_id FROM tcga.{table}
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
            geneStrings = [row[0] for row in rows]
            temp_geneLists = [genesString.split(',') for genesString in
                              geneStrings]
            geneLists.extend(temp_geneLists)
        return geneLists

    def _update_gene_coverage(self):
        """Populate object's gene_coverage dictionary."""
        total_patients = 0  # will hold total number of patients
        gene_coverage = dict()
        if self.yale_proj_ids:
            (total_temp, counts_temp) = self._get_gene_counts_yale()
            total_patients += total_temp
            for gene in counts_temp.iterkeys():
                if gene_coverage.has_key(gene):
                    gene_coverage[gene] += counts_temp[gene]
                else:
                    gene_coverage[gene] = counts_temp[gene]
        if self.tcga_proj_abbrvs:
            (total_temp, counts_temp) = self._get_gene_counts_tcga()
            total_patients += total_temp
            for gene in counts_temp.iterkeys():
                if gene_coverage.has_key(gene):
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

    def _get_gene_counts_yale(self):
        """Fetch dictionary: gene -> (int) number of patients with mutation. 
        Also appends patient names for each gene to geneMatrix object.
        Dictionary may be empty if no pathway genes were mutated."""
        total_patients = None
        count_dict = dict()
        projIdsTuple = self.yale_proj_ids
        cmd0 = """SET SESSION group_concat_max_len = 30000;"""
        cmd1 = """SELECT count(DISTINCT patient_id) AS num_patients FROM 
            (SELECT patient_id FROM mutations_tumor_normal NATURAL JOIN normals 
            NATURAL JOIN patients WHERE project_id IN ({projGroupStr}) 
            {patient_filter} AND `Variant_Classification` <> 'Silent'
            GROUP BY patient_id HAVING count(*) <= {max_mutations} ) p2;""". \
            format(max_mutations=self.max_mutations,
                   projGroupStr=','.join(str(i) for i in projIdsTuple),
                   patient_filter=self._build_patient_filter_str(form="AND"))
        # HUGO LIST AND PATIENT COUNTS
        cmd2 = """SELECT hugo_symbol, count(DISTINCT patient_id) AS n_patients, 
            GROUP_CONCAT(DISTINCT patient_id) AS patients
            FROM mutations_tumor_normal m NATURAL JOIN normals
            NATURAL JOIN patients
            # gene subset in pathway of interest
            INNER JOIN refs.`pathway_gene_link` pgl
                ON m.entrez_gene_id = pgl.entrez_id
            {expression_filter_pgl}
            NATURAL JOIN
            #good patients
            (SELECT patient_id FROM mutations_tumor_normal NATURAL JOIN normals 
                NATURAL JOIN patients WHERE project_id IN ({projGroupStr})
                {patient_filter} AND `Variant_Classification` <> 'Silent'
                GROUP BY patient_id HAVING count(*) <= {max_mutations} )  p
            WHERE path_id = {path_id} AND `Variant_Classification` <> 'Silent'
            {patient_filter}
            GROUP BY hugo_symbol
            ORDER BY hugo_symbol;""". \
            format(path_id=self.path_id,
                   max_mutations=self.max_mutations,
                   projGroupStr=','.join(
                       str(i) for i in projIdsTuple),
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
            total_patients = int(cur.fetchone()[0])
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
                raise Exception("Pathway coverage query gives inconsistent " +
                                "patient counts and patient names.")
            self.geneMatrix.addGenePatients(gene, patient_names)
            count_dict[gene] = coverage
        return total_patients, count_dict

    def _get_gene_counts_tcga(self):
        """ Fetch dictionary: gene -> int percentage of patients with mutation.
        Dictionary may be empty if no pathway genes were mutated."""
        total_patients = 0  # will increase as inspect tables
        count_dict = dict()
        rows = None
        for tcga_table in self.tcga_proj_abbrvs:
            cmd0 = """SET SESSION group_concat_max_len = 30000;"""
            cmd1 = """SELECT count(DISTINCT patient_id) AS n_patient FROM 
                (SELECT patient_id FROM tcga.{table} 
                WHERE `Variant_Classification` <> 'Silent' {patient_filter}
                GROUP BY patient_id HAVING count(*) <= {max_mutations}) p2;""" \
                .format(max_mutations=self.max_mutations, table=tcga_table,
                        patient_filter=self._build_patient_filter_str(
                            form="AND"))
            # HUGO LIST AND PATIENT COUNTS
            cmd2 = """SELECT hugo_symbol, count(DISTINCT patient_id)
                AS n_patients, GROUP_CONCAT(DISTINCT patient_id) AS patients
                FROM tcga.{table} t
                # gene subset in pathway of interest
                INNER JOIN refs.`pathway_gene_link` pgl
                ON t.entrez_id = pgl.entrez_id {expression_filter_pgl}
                NATURAL JOIN
                #good patients
                (SELECT patient_id FROM tcga.{table} t
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
                self.geneMatrix.addGenePatients(gene, patient_names)
                if count_dict.has_key(gene):
                    count_dict[gene] += coverage;
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
        for pway_size in xrange(1,
                                self.G):  # WAS xrange(1,self.G - self.max_mutations):
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

    def __init__(self, yale_proj_ids, tcga_proj_abbrvs, name_suffix=None):
        self.yale_proj_ids = yale_proj_ids
        self.tcga_proj_abbrvs = tcga_proj_abbrvs
        self.name_suffix = name_suffix
        self.root_name = self._get_root_filename(yale_proj_ids,
                                                 tcga_proj_abbrvs)

    def _get_root_filename(self, yale_proj_ids, tcga_proj_abbrvs):
        """Root file name for output files."""
        base_str = 'pathways_pvalues'
        yale_substr = str(yale_proj_ids).replace(',', '').replace(' ',
                                                                  '_').replace(
            '(', '').replace(')', '')
        if yale_proj_ids:
            yale_substr = '_' + yale_substr
        tcga_substr = '_'.join(tcga_proj_abbrvs)
        if tcga_proj_abbrvs:
            tcga_substr = '_' + tcga_substr
        root_name = base_str + yale_substr + tcga_substr
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

    def __init__(self, yale_proj_ids, tcga_proj_abbrvs, patient_ids=list(),
                 name_suffix=None, max_mutations=None, expressed_table=None,
                 ignore_genes=list()):
        # create self.root_name
        GenericPathwayFileProcessor.__init__(self, yale_proj_ids,
                                             tcga_proj_abbrvs,
                                             name_suffix=name_suffix)
        self.filter_patient_ids = patient_ids
        self.max_mutations = max_mutations
        self.expressed_table = expressed_table
        self.ignore_genes = ignore_genes

    def get_ordered_pway_list(self):
        file_name = self.root_name + '.txt'
        # max_lookup_rows = 100
        allPathways = list()
        with open(file_name, 'r') as file:
            for line in file:
                row = line.strip().split()
                path_id = int(row[0])
                pval = float(row[1])
                psize = int(row[2])
                peffect = int(row[3])
                runtime = float(row[4])
                # set up pathway object
                pway = PathwaySummary(path_id, self.yale_proj_ids,
                                      self.tcga_proj_abbrvs,
                                      patient_ids=self.filter_patient_ids,
                                      max_mutations=self.max_mutations,
                                      expressed_table=self.expressed_table,
                                      ignore_genes=self.ignore_genes)
                pway.set_up_from_file(pval, psize, peffect, runtime)
                allPathways.append(pway)
        # sort pathways by effect_size : actual_size
        allPathways.sort(key=lambda pway: self._get_size_ratio(pway.n_effective,
                                                               pway.n_actual),
                         reverse=True)
        # sort pathways by p-value
        allPathways.sort(key=lambda pway: pway.p_value, reverse=False)
        # # for pway in allPathways[0:max_lookup_rows+1]:
        # for pway in allPathways:
        # if pway.p_value < 0.1:
        # pway.update_exclusive_cooccurring_coverage()
        return allPathways

    def _get_size_ratio(self, n_effective, n_actual):
        """Ratio used for initial sorting. >1 if large effective.
        0<r<1 if small_effective. 0 if n_actual is 0."""
        if n_actual:
            return float(n_effective) / n_actual
        else:
            return 0


class PathwayDetailedFileWriter(GenericPathwayFileProcessor):
    """Writes detailed postprocessing file:
    pathway names, pvalues and gene info."""

    def __init__(self, yale_proj_ids, tcga_proj_abbrvs, pway_object_list,
                 name_suffix=None):
        # create self.root_name
        GenericPathwayFileProcessor.__init__(self, yale_proj_ids,
                                             tcga_proj_abbrvs,
                                             name_suffix=name_suffix)
        self.allPathways = pway_object_list
        self.nameDict = self.getPathwayNameDict()
        self.outfile_name = self.root_name + '_pretty.txt'
        self.matrix_folder = 'matrix_txt'

    def getPathwayNameDict(self):
        """Gets name for all pathways, stored in dict: pathid -> pathname."""
        rows = None
        pathwayDict = dict()
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
                    % (row_count))
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
            pathwayDict[path_id] = path_name
        return pathwayDict

    def dictToStructure(self, coverage_dict):
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
                coverage_string = self.dictToStructure(pway.gene_coverage)
                pway.exclusive_genes.sort(key=lambda
                    gene: pway.gene_coverage[gene], reverse=True)
                pway.cooccurring_genes.sort(key=lambda
                    gene: pway.gene_coverage[gene], reverse=True)
                exclusive_string = '{' + ','.join(
                    [repr(i) for i in pway.exclusive_genes]) + '}'
                cooccurring_string = '{' + ','.join(
                    [repr(i) for i in pway.cooccurring_genes]) + '}'
                out.write(
                    "{path_id}\t{name}\t{n_actual}\t{n_effective}\t{p_value:.3e}\t" +
                    "{runtime:.2f}\t{exclusive_string}\t{cooccurring_string}\t" +
                    "{coverage_string}\n"
                    .format(path_id=pway.path_id,
                            name=path_name,
                            coverage_string=coverage_string,
                            n_actual=pway.n_actual,
                            n_effective=pway.n_effective,
                            p_value=pway.p_value,
                            runtime=pway.runtime,
                            exclusive_string=exclusive_string,
                            cooccurring_string=cooccurring_string))

    def write_matrix_files(self, pway):
        """Write text file containing presence matrix for patient-gene pair."""
        if not os.path.exists(self.matrix_folder):
            os.mkdir(self.matrix_folder)
        matrix_filename = self.matrix_folder + os.sep + self.root_name + '_matrix' + str(
            pway.path_id) + '.txt'
        with open(matrix_filename, 'w') as outfile:
            pway.geneMatrix.export_matrix(outfile, pway.exclusive_genes)


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

    def _fetch_expressed_genome_size(self, expressed_table):
        """Count genes via SQL query: assumes row count equals gene count."""
        cmd1 = """SELECT count(*) FROM tcga.{table_name};""".format(
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


def main():
    """Arguments: projIds --patients patient_file."""

    parser = argparse.ArgumentParser()
    parser.add_argument("proj_ids", help="list of project identifiers",
                        nargs="+")  # ids are strings, e.g. ["1","22","LUSC"]
    parser.add_argument("-p", "--patients", type=file, dest="patients_file",
                        help="file containing patients to include")
    parser.add_argument("-g", "--genes", nargs="+",
                        help="ignore pathways that don't contain these genes " +
                             "(given by hugo symbol)")
    parser.add_argument("-gs", "--genome",
                        help="limit genome size to protein-coding genes.",
                        choices=['no_pseudo', 'anything', 'protein-coding',
                                 'inc_misc_chr'])
    parser.add_argument("-s", "--suffix",
                        help="optional descriptive string to append to filenames.")
    parser.add_argument("--cutoff", type=int, default=500,
                        help="maximum number of mutations allowed for sample inclusion.")
    parser.add_argument("--expression", nargs=1,
                        help="name of table containing entrez_id for expressed genes.")
    parser.add_argument("--ignore", nargs="+",
                        help="list of genes (by hugo symbol) to ignore when calculating p-values (and looking up pathway sizes).")
    args = parser.parse_args()

    # max_mutations
    max_mutations = args.cutoff

    if args.expression:
        args.expression = args.expression[0]

    ignore_genes = args.ignore

    genome_size = BackgroundGenomeFetcher(args.genome,
                                          args.expression).genome_size

    # get patient list. maybe empty list.
    patient_list = list()
    if args.patients_file:
        for line in args.patients_file:
            temp_line = line.strip('\n')
            if temp_line.isdigit():
                patient_list.append(int(temp_line))  # integer if possible
            else:
                patient_list.append(temp_line)  # strings for tcga projects
        print("Loaded {} patients.".format(len(patient_list)))
    else:
        print('No patients file provided. Using all patients.')

    # split projIds into yale, tcga    
    yale_proj_ids = tuple(int(i) for i in args.proj_ids if i.isdigit())
    tcga_proj_abbrvs = tuple(i for i in args.proj_ids if not i.isdigit())
    print("Yale projects: {}".format(str(yale_proj_ids)))
    print("TCGA projects: {}".format(str(tcga_proj_abbrvs)))

    # get genes of interest, if any
    if args.genes:
        interest_genes = tuple(args.genes)
    else:
        interest_genes = tuple()
    print("Interest genes: {}".format(str(interest_genes)))

    idFetcher = PathwayIdsFetcher(interest_genes)
    all_path_ids = idFetcher.fetchPathwayIds()

    for pathway_number in all_path_ids:
        # Populate pathway object, and time pvalue calculation
        start = timeit.default_timer()
        pway = PathwaySummary(pathway_number, yale_proj_ids, tcga_proj_abbrvs,
                              patient_ids=patient_list,
                              max_mutations=max_mutations,
                              expressed_table=args.expression,
                              ignore_genes=ignore_genes)
        pway.set_pathway_size()
        pway.populate_patient_info()
        lcalc = LCalculator(pway, genome_size)  # include optional genome_size
        lcalc.run()
        runtime = timeit.default_timer() - start
        # Write results to 'basic' file
        basicWriter = PathwayBasicFileWriter(yale_proj_ids, tcga_proj_abbrvs,
                                             name_suffix=args.suffix)
        basicWriter.write_pvalue_file(lcalc, runtime)

    # Gather all pathway stats from text file
    assembler = PathwayListAssembler(yale_proj_ids, tcga_proj_abbrvs,
                                     patient_ids=patient_list,
                                     name_suffix=args.suffix,
                                     max_mutations=max_mutations,
                                     expressed_table=args.expression,
                                     ignore_genes=ignore_genes)
    pway_list = assembler.get_ordered_pway_list()

    # Rank pathways, gather extra stats and write to final file
    final_writer = PathwayDetailedFileWriter(yale_proj_ids,
                                             tcga_proj_abbrvs, pway_list,
                                             name_suffix=args.suffix)
    final_writer.write_detailed_file()


if __name__ == '__main__':
    dbvars = {'host': 'localhost', 'db': 'CancerDB',
              'read_default_file': "~/.my.cnf"}
    main()

    

