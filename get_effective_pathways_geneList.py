#!/Users/Stephen/Library/Enthought/Canopy_64bit/User/bin/python

import MySQLdb as mdb
from scipy.misc import comb
from numpy import *
from scipy import stats
import timeit
import sys
import argparse
# import pdb

class NonSingleResult(Exception):
    pass


class Patient():
    def __init__(self, patient_id, n_mutated, is_mutated, is_yale):
        self.patient_id = patient_id
        self.n_mutated = n_mutated
        self.is_mutated = is_mutated
        self.is_yale = is_yale


class PathwaySummary():
    """Holds pathway information, and can fetch info from db."""
    def __init__(self,pathway_number, yale_proj_ids, tcga_proj_abbrvs):
        self.path_id = pathway_number
        self.n_actual = None
        self.patients = list() # tuples. (patient_id, n_mutations, is_mutated)
        self.n_effective = None
        self.p_value = None
        self.max_mutations = 500
        self.yale_proj_ids = yale_proj_ids
        self.tcga_proj_abbrvs = tcga_proj_abbrvs
        # populated during post-processing:
        self.gene_coverage = dict()
        self.exclusive_genes = list()
        self.cooccurring_genes = list()
        self.runtime = None  #only set up by file reader

    def set_up_from_file(self,pval,psize,peffect,runtime):
        """Manually specify pathway summary properties."""
        self.p_value = pval
        self.n_actual = psize
        self.n_effective = peffect
        self.runtime = runtime
        # # fetch gene_coverage, exclusive_genes, cooccurring genes
        # self._populate_exclusive_cooccurring()
        # self._update_gene_coverage()

    def set_pathway_size(self):
        """Lookup pathway size in DB. Save to data attribute."""
        pway_size = None

        # GET pway_size
        cmd1 = """SELECT count(DISTINCT entrez_id) AS pway_size 
        FROM refs.pathway_gene_link WHERE entrez_id IS NOT NULL 
        AND path_id = {pathid} GROUP BY path_id;""".format(
            pathid=self.path_id)  
        
        try:
            con = mdb.connect(**dbvars)
            cur = con.cursor()
            cur.execute(cmd1)
            rowCount = cur.rowcount
            if not rowCount == 1:
                raise NonSingleResult("Result contains %g rows Ids for pathway %s." 
                    % (rowCount, pathway_number))
            pway_size = cur.fetchone()
        except mdb.Error as e:
            print "Error %d: %s" % (e.args[0],e.args[1])
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

    def _append_yale_patients(self):
        """Get patient-pathway gene overlap info from databse.
        # Returns tuple collection of row tuples.
        Row tuple: (PATIENT_ID, N_PATIENT, BOOL_MUTATED), e.g. (678L, 323L, 1L)
        """
        rows = None
        
        ## GET mutated boolean for genes in specified pathway for GOOD 
        #  patients (<max_mutations), even those without somatic mutation. 
        cmd2 = """SELECT p.patient_id, n_patient, g.patient_id IS NOT NULL AS mutated FROM
            # good patients for project, mutation_count
            (SELECT patient_id, count(DISTINCT entrez_gene_id) AS n_patient FROM 
                mutations_tumor_normal NATURAL JOIN normals NATURAL JOIN patients 
            WHERE project_id IN ({projGroupStr})
            GROUP BY patient_id HAVING count(*) <= {max_mutations})  p 
        LEFT JOIN
            # patients in project with mutation in pathway
            (SELECT patient_id FROM 
                    (SELECT entrez_id FROM refs.`pathway_gene_link` pwg 
                        WHERE path_id={path_id}) pwg 
                INNER JOIN mutations_tumor_normal m ON pwg.entrez_id=m.entrez_gene_id 
                NATURAL JOIN normals NATURAL JOIN patients 
                WHERE project_id IN ({projGroupStr}) 
                GROUP BY patient_id) g 
        ON p.patient_id = g.patient_id;""".format(path_id=self.path_id, 
            max_mutations=self.max_mutations, 
            projGroupStr = ','.join(str(i) for i in self.yale_proj_ids))

        
        try:
            con = mdb.connect(**dbvars)
            cur = con.cursor()
            cur.execute(cmd2)
            rows = cur.fetchall()
        except mdb.Error as e:
            print "Error %d: %s" % (e.args[0],e.args[1])
        finally:    
            if con: 
                con.close()

        for row in rows:
            patient_id = row[0]
            n_mutated = row[1]
            is_mutated = row[2]>0
            is_yale = True
            patient = Patient(patient_id, n_mutated, is_mutated, is_yale)
            self.patients.append(patient)

    def _append_tcga_patients(self):
        """Gets patient-pathway gene overlap info from databse.
        # Returns tuple collection of row tuples.
        Row tuple: (PATIENT_ID, N_PATIENT, BOOL_MUTATED), e.g. (678L, 323L, 1L)
        """
        for tcga_table in self.tcga_proj_abbrvs:
            rows = None
            
            ## GET mutated boolean for genes in specified pathway for GOOD 
            #  patients (<max_mutations), even those without somatic mutation. 
            cmd2 = """SELECT p.patient_id, n_patient, g.patient_id IS NOT NULL AS mutated FROM
                # good patients, mutation_count
                (SELECT patient_id, count(DISTINCT entrez_id) AS n_patient FROM tcga.{table} 
                    GROUP BY patient_id HAVING count(*) <= {max_mutations}) p
                LEFT JOIN
                # pathway mutation counts above 0
                (SELECT patient_id FROM
                    # get distinct genes from patients
                    (SELECT DISTINCT patient_id, entrez_id FROM tcga.{table}) pg
                    # join to pathway genes
                        INNER JOIN 
                    (SELECT entrez_id FROM refs.pathway_gene_link WHERE path_id = {path_id}) pwg 
                    ON pwg.entrez_id = pg.entrez_id 
                    GROUP BY `patient_id`) g
                ON p.patient_id = g.patient_id;""".format(path_id=self.path_id, 
                max_mutations=self.max_mutations,table=tcga_table)
            
            try:
                con = mdb.connect(**dbvars)
                cur = con.cursor()
                cur.execute(cmd2)
                rows = cur.fetchall()
            
            except mdb.Error as e:
                print "Error %d: %s" % (e.args[0],e.args[1])
            finally:    
                if con: 
                    con.close()
            
            for row in rows:
                patient_id = row[0]
                n_mutated = row[1]
                is_mutated = row[2]>0
                is_yale = False
                patient = Patient(patient_id, n_mutated, is_mutated, is_yale)
                self.patients.append(patient)

    def update_exclusive_cooccurring_coverage(self):
        """ POSTPROCESSING. Gather pathway gene info and write detailed output."""
        self._populate_exclusive_cooccurring()
        self._update_gene_coverage()

    def _populate_exclusive_cooccurring(self):
        """ Postprocessing step. Look up gene combinations hit, sort into sets."""
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
            if len(geneList)==1:
                exclusive_gene_set.add(geneList[0])
            for gene in geneList:
                all_hit_genes_set.add(gene)
        cooccurring_gene_set = all_hit_genes_set.difference(exclusive_gene_set)

        onlyPairedWithExclusive = dict(zip(cooccurring_gene_set,
            [True for i in cooccurring_gene_set]))
        for geneList in gene_combs_list:
            if not exclusive_gene_set.intersection(set(geneList)):
                for gene in geneList:
                    onlyPairedWithExclusive[gene] = False

        tagAlongGenes = [gene for gene in cooccurring_gene_set 
            if onlyPairedWithExclusive[gene]]
        coExclusiveGenes = [gene for gene in cooccurring_gene_set if not 
            onlyPairedWithExclusive[gene]]

        for gene in coExclusiveGenes:
            exclusive_gene_set.add(gene)

        self.exclusive_genes = sorted(list(exclusive_gene_set))
        self.cooccurring_genes = sorted(tagAlongGenes)
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
            SELECT patient_id, group_concat(DISTINCT hugo_symbol ORDER BY hugo_symbol
                 SEPARATOR ',') AS symbols 
            FROM 
                # gene subset in pathway of interest            
                (SELECT entrez_id FROM refs.`pathway_gene_link` WHERE path_id = {path_id}) pgl
            INNER JOIN mutations_tumor_normal m 
            ON m.entrez_gene_id = pgl.entrez_id
            NATURAL JOIN normals 
            NATURAL JOIN patients p
            WHERE project_id IN ({projGroupStr})
            GROUP BY patient_id
        ) g
        INNER JOIN
            #good patients
            (SELECT patient_id FROM mutations_tumor_normal NATURAL JOIN normals 
                NATURAL JOIN patients WHERE project_id IN ({projGroupStr})
                GROUP BY patient_id HAVING count(*) <= {max_mutations} )  p2
        ON g.patient_id = p2.patient_id;""".format(path_id=self.path_id,  
            max_mutations=self.max_mutations,
            projGroupStr = ','.join(str(i) for i in projIdsTuple))

        try:
            con = mdb.connect(**dbvars)
            cur = con.cursor()
            cur.execute(cmd)
            rows = cur.fetchall()
        except mdb.Error as e:
            print "Error %d: %s" % (e.args[0],e.args[1])
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
                SELECT patient_id, group_concat(DISTINCT hugo_symbol ORDER BY hugo_symbol
                     SEPARATOR ',') AS symbols
                FROM tcga.{table} t
                # gene subset in pathway of interest
                INNER JOIN refs.`pathway_gene_link` pgl ON t.entrez_id = pgl.entrez_id
                NATURAL JOIN
                #good patients
                (SELECT patient_id, count(DISTINCT entrez_id) AS n_patient 
                    FROM tcga.{table} GROUP BY patient_id HAVING count(*) <= {max_mutations}) p
                WHERE path_id = {path_id}
                GROUP BY patient_id
            ) g;""".format(path_id=self.path_id, max_mutations=self.max_mutations,
            table=tcga_table)

            try:
                con = mdb.connect(**dbvars)
                cur = con.cursor()
                cur.execute(cmd)
                rows = cur.fetchall()
            except mdb.Error as e:
                print "Error %d: %s" % (e.args[0],e.args[1])
            finally:    
                if con: 
                    con.close()

            geneStrings = [row[0] for row in rows]
            temp_geneLists = [genesString.split(',') for genesString in geneStrings]
            geneLists.extend(temp_geneLists)
        return geneLists

    def _update_gene_coverage(self):
        """ Populate object's gene_coverage dictionary."""

        t_patients = 0  # will hold total number of patients
        gene_coverage = dict()

        if self.yale_proj_ids:
            (t_temp, counts_temp) = self._get_gene_counts_yale()
            t_patients += t_temp
            for gene in counts_temp.iterkeys():
                if gene_coverage.has_key(gene):
                    gene_coverage[gene] += counts_temp[gene]
                else:
                    gene_coverage[gene] = counts_temp[gene]

        if self.tcga_proj_abbrvs:
            (t_temp, counts_temp) = self._get_gene_counts_tcga()
            t_patients += t_temp
            for gene in counts_temp.iterkeys():
                if gene_coverage.has_key(gene):
                    gene_coverage[gene] += counts_temp[gene]
                else:
                    gene_coverage[gene] = counts_temp[gene]

        # convert to coverage by dividing by number of patients
        for gene in gene_coverage:
            gene_coverage[gene] = float(gene_coverage[gene])/t_patients*100
        
        self.gene_coverage = gene_coverage
        return

    def _get_gene_counts_yale(self):
        """ Fetch dictionary: gene -> (int) number of patients with mutation. 
        Dictionary may be empty if no pathway genes were mutated."""
        t_patients = None
        count_dict = dict()
        projIdsTuple = self.yale_proj_ids
        
        cmd1 = """SELECT count(DISTINCT patient_id) AS num_patients FROM 
            (SELECT patient_id FROM mutations_tumor_normal NATURAL JOIN normals 
            NATURAL JOIN patients WHERE project_id IN ({projGroupStr}) 
            GROUP BY patient_id HAVING count(*) <= {max_mutations} ) p2;""".format(
            max_mutations=self.max_mutations,
            projGroupStr = ','.join(str(i) for i in projIdsTuple))

        # HUGO LIST AND PATIENT COUNTS
        cmd2 = """SELECT hugo_symbol, count(DISTINCT patient_id) AS n_patients
            FROM mutations_tumor_normal m NATURAL JOIN normals NATURAL JOIN patients
            # gene subset in pathway of interest
            INNER JOIN refs.`pathway_gene_link` pgl ON m.entrez_gene_id = pgl.entrez_id
            NATURAL JOIN
            #good patients
            (SELECT patient_id FROM mutations_tumor_normal NATURAL JOIN normals 
                NATURAL JOIN patients WHERE project_id IN ({projGroupStr})
                GROUP BY patient_id HAVING count(*) <= {max_mutations} )  p
            WHERE path_id = {path_id}
            GROUP BY hugo_symbol
            ORDER BY hugo_symbol;""".format(path_id=self.path_id,
                max_mutations=self.max_mutations,
                projGroupStr = ','.join(str(i) for i in projIdsTuple))
        
        try:
            con = mdb.connect(**dbvars)
            cur = con.cursor()
            
            cur.execute(cmd1)
            rowCount = cur.rowcount
            if not rowCount == 1:
                raise NonSingleResult("Result contains %g rows Ids for pathway %s."
                     % (rowCount, pathway_number))
            
            t_patients = int(cur.fetchone()[0])

            cur.execute(cmd2)
            rowCount = cur.rowcount
            if not rowCount:
                # NO GENES MUTATED. n_effective < n_pathway
                return t_patients, count_dict
            
            rows = cur.fetchall()
        
        except mdb.Error as e:
            print "Error %d: %s" % (e.args[0],e.args[1])
        finally:    
            if con: 
                con.close()
        
        for geneRow in rows:
            gene = geneRow[0]
            coverage = int(geneRow[1])
            count_dict[gene] = coverage

        return t_patients, count_dict

    def _get_gene_counts_tcga(self):
        """ Fetch dictionary: gene -> (int) percentage of patients with mutation. 
        Dictionary may be empty if no pathway genes were mutated."""
        t_patients = 0  # will increase as inspect tables
        count_dict = dict()
        rows = None
        
        for tcga_table in self.tcga_proj_abbrvs:

            cmd1 = """SELECT count(DISTINCT patient_id) AS n_patient FROM 
                (SELECT patient_id FROM tcga.{table} 
                GROUP BY patient_id HAVING count(*) <= {max_mutations}) p2;""".format(
                max_mutations=self.max_mutations,table=tcga_table)

            # HUGO LIST AND PATIENT COUNTS
            cmd2 = """SELECT hugo_symbol, count(DISTINCT patient_id) AS n_patients
                FROM tcga.{table} t
                # gene subset in pathway of interest
                INNER JOIN refs.`pathway_gene_link` pgl ON t.entrez_id = pgl.entrez_id
                NATURAL JOIN
                #good patients
                (SELECT patient_id FROM tcga.{table} t
                    GROUP BY patient_id HAVING count(*) <= {max_mutations} )  p
                WHERE path_id = {path_id}
                GROUP BY hugo_symbol
                ORDER BY hugo_symbol;""".format(path_id=self.path_id,
                    max_mutations=self.max_mutations,table=tcga_table)
            
            try:
                con = mdb.connect(**dbvars)
                cur = con.cursor()
                
                cur.execute(cmd1)
                rowCount = cur.rowcount
                if not rowCount == 1:
                    raise NonSingleResult("Result contains %g rows Ids for pathway %s."
                         % (rowCount, pathway_number))
                
                t_patients += int(cur.fetchone()[0])

                cur.execute(cmd2)
                rowCount = cur.rowcount
                if not rowCount:
                    # NO GENES MUTATED. n_effective < n_pathway
                    return t_patients, count_dict
                
                rows = cur.fetchall()
            
            except mdb.Error as e:
                print "Error %d: %s" % (e.args[0],e.args[1])
            finally:    
                if con: 
                    con.close()
            
            for geneRow in rows:
                gene = geneRow[0]
                coverage = int(geneRow[1])
                if count_dict.has_key(gene):
                    count_dict[gene] += coverage; 
                else:    
                    count_dict[gene] = coverage

        return t_patients, count_dict


class LCalculator():
    """Calculates likelihood of observing pathway mutations in patients, and 
    MLE pathway size from these observations. Takes a pathwaySummary object."""

    def __init__(self, pway):

        self.G = 18852 # genes in genome
        self.pway = pway
        self.likelihood = None
        self.ne = None
        self.ne_ll = None
        self.D = None
        self.pvalue = None

    def run(self):
        """ Calculate likelihood and maximum likelihood estimate."""
        self.likelihood = self._get_pway_likelihood(self.pway.n_actual)
        (ne,lastll) = self._get_ne()
        self.ne = ne
        self.ne_ll = lastll
        self.D = -2*self.likelihood + 2*self.ne_ll
        self.pvalue = 1-stats.chi2.cdf(self.D,1)

        self.pway.n_effective = self.ne
        self.pway.p_value = self.pvalue

    def _get_pway_likelihood(self, pway_size=None):
        if pway_size is None:
            pway_size = self.pway_size
            
        prob_list = list()
        #iterate over patients
        for patient in self.pway.patients:
            n_patient = patient.n_mutated
            mutated = patient.is_mutated
            p_no_mut = exp(math.log(comb(self.G - pway_size, n_patient, 
                exact=True)) - math.log(comb(self.G, n_patient, exact=True)))
            if mutated:
                p = 1-p_no_mut
            else:
                p = p_no_mut
            prob_list.append(log(p))
        prob_array = array(prob_list)
        return prob_array.sum()
        
    def _get_ne(self):
        last_ll = None
        #improved = False
        ne = None
        #profile = list()

        #check last 2 vals to check for decline:
        penult = self._get_pway_likelihood(pway_size=self.G - 500 - 1)
        ult = self._get_pway_likelihood(pway_size=self.G - 500)
        if ult > penult:
            ne = self.G
            return (ne,ult)

        # at this stage, there will be a max before Genome size
        for pway_size in range(1,self.G - 500):
            this_ll = self._get_pway_likelihood(pway_size=pway_size)
            #profile.append(this_ll)
            #if mod(pway_size,100)==0:
            #print this_ll

            if last_ll is None:
                last_ll = this_ll
                
            #if improved is False and (this_ll > self.likelihood or 
            #    pway_size >= self.pway_size):
            #    improved = True
                
            if this_ll < last_ll:
                ne = pway_size - 1
                
                break
            
            last_ll = this_ll
        return (ne,last_ll)


class GenericPathwayFileProcessor():
    """Generic object that can convert yale_proj_ids and tcga_proj_abbrvs 
    to file_name."""
    def __init__(self, yale_proj_ids, tcga_proj_abbrvs):
        self.root_name = self._get_root_filename(yale_proj_ids, tcga_proj_abbrvs)
        self.yale_proj_ids = yale_proj_ids
        self.tcga_proj_abbrvs = tcga_proj_abbrvs
        
    def _get_root_filename(self,yale_proj_ids, tcga_proj_abbrvs):
        """Root file name for output files."""
        
        base_str = 'pathways_pvalues'
        yale_substr = str(yale_proj_ids).replace(',','').replace(' ','_').replace(
            '(','').replace(')','')
        if yale_proj_ids:
            yale_substr = '_' + yale_substr
        tcga_substr = '_'.join(tcga_proj_abbrvs)
        if tcga_proj_abbrvs:
            tcga_substr = '_' + tcga_substr
        root_name = base_str + yale_substr + tcga_substr
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
                path_id,lcalc.pvalue,lcalc.pway.n_actual,lcalc.ne,runtime))


class PathwayListAssembler(GenericPathwayFileProcessor):
    """Builds ordered list of pathways from basic p-value file."""

    def get_ordered_pway_list(self):
        file_name = self.root_name + '.txt'
        #max_lookup_rows = 100

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
                    self.tcga_proj_abbrvs)
                pway.set_up_from_file(pval,psize,peffect,runtime)

                allPathways.append(pway)
        
        # sort pathways by effect_size : actual_size
        allPathways.sort(key=lambda pway: float(pway.n_effective)/pway.n_actual,
            reverse=True)
        # sort pathways by p-value
        allPathways.sort(key=lambda pway: pway.p_value,reverse=False)

        ## for pway in allPathways[0:max_lookup_rows+1]:
        # for pway in allPathways:
        #     if pway.p_value < 0.1:
        #         pway.update_exclusive_cooccurring_coverage()

        return allPathways


class PathwayDetailedFileWriter(GenericPathwayFileProcessor):
    """Writes detailed postprocessing file: pathway names, pvalues and gene info."""
    def __init__(self, yale_proj_ids, tcga_proj_abbrvs, pway_object_list):
        # create self.root_name
        GenericPathwayFileProcessor.__init__(self,yale_proj_ids, tcga_proj_abbrvs)
        self.allPathways = pway_object_list
        self.nameDict = self.getPathwayNameDict()
        self.outfile_name = self.root_name + '_pretty.txt'
        
    def getPathwayNameDict(self):
        """Gets name for all pathways, stored in dictionary: pathid -> pathname."""
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
            rowCount = cur.rowcount
            if not rowCount > 1:
                raise Exception("Result contains %g rows Ids for pathway lookup." 
                    % (rowCount))
            
            rows = cur.fetchall()
        
        except mdb.Error as e:
            print "Error %d: %s" % (e.args[0],e.args[1])
        finally:    
            if con: 
                con.close()
        
        # rows is [[id,name],[id,name],...]
        for pair in rows:
            path_id = int(pair[0])
            path_name = pair[1]
            pathwayDict[path_id] = path_name

        return pathwayDict

    def dictToStructure(self,coverage_dict):
        """ Get matlab command to convert dictionary to structure.
        e.g. Dict: 'BRAF' -> 34. --> struct('BRAF',34)."""
        
        # struct_string = repr(coverage_dict) # "('BRAF':34, 'KRAS':16)"
        # struct_string = struct_string.replace(":",",")
        # struct_string = struct_string.replace("{","(").replace("}",")")
        # struct_string = "struct" + struct_string
        # return struct_string
        
        # name,coverage pairs
        pairs = ",".join(["{!r},{:.2f}".format(gene,coverage_dict[gene]) for 
            gene in coverage_dict])
        struct_string = "struct(" + pairs + ")"
        return struct_string





    def write_detailed_file(self):
        """Perform lookup of coverage etc for low pvalue pathways and write 
        ordered list of pathways plus info to file."""       
        bufsize = 1
        with open(self.outfile_name,'w', bufsize) as out:
            for pway in self.allPathways:
                # Get extra info, if p_value is low
                if pway.p_value < 0.1:
                    pway.update_exclusive_cooccurring_coverage()

                path_name = self.nameDict[pway.path_id]
                coverage_string = self.dictToStructure(pway.gene_coverage)

                pway.exclusive_genes.sort(key=lambda 
                    gene: pway.gene_coverage[gene], reverse=True)
                pway.cooccurring_genes.sort(key=lambda 
                    gene: pway.gene_coverage[gene],reverse=True)

                exclusive_string = '{' + ','.join(
                    [repr(i) for i in pway.exclusive_genes]) + '}'
                cooccurring_string = '{' + ','.join(
                    [repr(i) for i in pway.cooccurring_genes]) + '}'

                out.write(("{path_id}\t{name}\t{n_actual}\t{n_effective}\t{p_value:.3e}\t"+
                    "{runtime:.2f}\t{exclusive_string}\t{cooccurring_string}\t"+
                    "{coverage_string}\n").format(path_id=pway.path_id, name=path_name,
                    coverage_string=coverage_string,n_actual=pway.n_actual,
                    n_effective=pway.n_effective,p_value=pway.p_value,
                    runtime=pway.runtime, exclusive_string=exclusive_string,
                    cooccurring_string=cooccurring_string))        



class PathwayIdsFetcher():
    def __init__(self, interest_genes):
        self.interest_genes = interest_genes

    def fetchPathwayIds(self):
        """Get pathway ids containing genes in (possibly empty) interest set."""
        rows = None
        all_path_ids = list()

        if self.interest_genes:
            genes_string = repr(tuple(self.interest_genes))
            genes_string = genes_string.replace(",)",")")
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
                raise Exception("Result contains %g rows Ids for pathway lookup." 
                    % (rowCount))
            
            rows = cur.fetchall()
        
        except mdb.Error as e:
            print "Error %d: %s" % (e.args[0],e.args[1])
        finally:    
            if con: 
                con.close()
        
        # rows is [[id,name],[id,name],...]
        for id in rows:
            all_path_ids.append(int(id[0]))
            
        return all_path_ids


def main():        
    """Calculate all effective pathway sizes and p-values for project(s) 
    provided as list of tuples."""
    proj_ids = eval(sys.argv[1])  # e.g. "(1,)" or "(1,2,'LUAD')"
    interest_genes = eval(sys.argv[2])  # e.g. "('CBL',)"

    # proj_ids = (22,'LUAD')
    # interest_genes = ('CBL',)

    yale_proj_ids = list()
    tcga_proj_abbrvs = list() #empty list of no tcga arguments

    for id in proj_ids:
        if isinstance(id,int):
            yale_proj_ids.append(id)
        elif isinstance(id,str):
            tcga_proj_abbrvs.append(id)
        else:
            raise Exception("Unrecognised project id in tuple.")

    yale_proj_ids = tuple(yale_proj_ids)
    tcga_proj_abbrvs = tuple(tcga_proj_abbrvs)

    #all_path_ids = range(1,1321)
    idFetcher = PathwayIdsFetcher(interest_genes)

    all_path_ids = idFetcher.fetchPathwayIds()

    #os.chdir('/Users/Stephen/Dropbox/Townsend/pathway_enrichment/')
    # for projIds in groups:

    for pathway_number in all_path_ids:

        start = timeit.default_timer()
        
        pway = PathwaySummary(pathway_number,yale_proj_ids,tcga_proj_abbrvs)
        pway.set_pathway_size()
        pway.populate_patient_info()
        lcalc = LCalculator(pway)
        lcalc.run()
        
        runtime = timeit.default_timer() - start
        basicWriter = PathwayBasicFileWriter(yale_proj_ids,tcga_proj_abbrvs)
        basicWriter.write_pvalue_file(lcalc, runtime)

    # Gather all pathway stats from text file
    assembler = PathwayListAssembler(yale_proj_ids, tcga_proj_abbrvs)
    pway_list = assembler.get_ordered_pway_list()

    # Rank pathways, gather extra stats and write to final file
    final_writer = PathwayDetailedFileWriter(yale_proj_ids, 
        tcga_proj_abbrvs, pway_list)
    final_writer.write_detailed_file()

def main2():        
    """Calculate all effective pathway sizes and p-values for project(s) 
    provided as list of tuples."""

    yale_proj_ids = (23,)
    tcga_proj_abbrvs = tuple()

    #os.chdir('/Users/Stephen/Dropbox/Townsend/pathway_enrichment/')

    # Gather all pathway stats from text file
    assembler = PathwayListAssembler(yale_proj_ids, tcga_proj_abbrvs)
    pway_list = assembler.get_ordered_pway_list()

    # Rank pathways, gather extra stats and write to final file
    final_writer = PathwayDetailedFileWriter(yale_proj_ids, 
        tcga_proj_abbrvs, pway_list)
    final_writer.write_detailed_file()

def runAll23():

    yale_proj_ids = (23,)
    tcga_proj_abbrvs = tuple()
    interest_genes = tuple()

    idFetcher = PathwayIdsFetcher(interest_genes)

    all_path_ids = idFetcher.fetchPathwayIds()

    #os.chdir('/Users/Stephen/Dropbox/Townsend/pathway_enrichment/')
    # for projIds in groups:

    for pathway_number in all_path_ids:

        start = timeit.default_timer()
        
        pway = PathwaySummary(pathway_number,yale_proj_ids,tcga_proj_abbrvs)
        pway.set_pathway_size()
        pway.populate_patient_info()
        lcalc = LCalculator(pway)
        lcalc.run()
        
        runtime = timeit.default_timer() - start
        basicWriter = PathwayBasicFileWriter(yale_proj_ids,tcga_proj_abbrvs)
        basicWriter.write_pvalue_file(lcalc, runtime)

    # Gather all pathway stats from text file
    assembler = PathwayListAssembler(yale_proj_ids, tcga_proj_abbrvs)
    pway_list = assembler.get_ordered_pway_list()

    # Rank pathways, gather extra stats and write to final file
    final_writer = PathwayDetailedFileWriter(yale_proj_ids, 
        tcga_proj_abbrvs, pway_list)
    final_writer.write_detailed_file()


def run_with_patient_subset():
    """Arguments: projIds --patients patient_file."""

    parser = argparse.ArgumentParser()
    parser.add_argument("proj_ids", help="list of project identifiers", 
        nargs="+") # ids are strings, e.g. ["1","22","LUSC"]
    parser.add_argument("-p", "--patients", type=file, dest="patients_file",
                        help="file containing patients to include")
    parser.add_argument("-g", "--genes", nargs="+",
                        help="file containing patients to include")
    args = parser.parse_args()
    
    # get patient list. maybe empty list.
    patient_list = list()
    if args.patients_file:
        print('Using patients file.')
        for line in patients_file:
            patient_list.append(line.strip('\n'))
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

    # idFetcher = PathwayIdsFetcher(interest_genes)
    # all_path_ids = idFetcher.fetchPathwayIds()

    



if __name__ == '__main__':
    dbvars = {'host':'localhost','db':'CancerDB','read_default_file':"~/.my.cnf"}
    # main()
    # main2()    
    # runAll23()
    run_with_patient_subset()

    


