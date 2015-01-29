__author__ = 'sgg'
from app import dbvars
import MySQLdb as mdb
from collections import defaultdict


def lookup_path_sizes_global():
    """Get pathway sizes when exclude genes are provided."""
    cmd = """SELECT path_id, count(DISTINCT entrez_id)
    FROM refs.pathway_gene_link pgl
    INNER JOIN refs.ncbi_entrez n ON pgl.entrez_id = n.geneId
    GROUP BY path_id;""".format()
    size_dict = dict()
    try:
        con = mdb.connect(**dbvars)
        cur = con.cursor()
        cur.execute(cmd)
        # assert isinstance(cur.rowcount, int)
        row_count = cur.rowcount
        if not row_count:
            print "No pathways found."
            return size_dict
        for row_no in xrange(row_count):
            row = cur.fetchone()
            size_dict[row[0]] = row[1]
    except mdb.Error as e:
        print "Error %d: %s" % (e.args[0], e.args[1])
    finally:
        if con:
            con.close()
    return size_dict


def lookup_path_sizes_exclude(ignore_genes):
    """Get pathway sizes when exclude genes are provided."""

    genes_string = repr(tuple(ignore_genes)).replace(",)", ")")
    genes_string = "WHERE symbol NOT IN {}".format(genes_string)

    cmd = """SELECT path_id, count(DISTINCT entrez_id)
    FROM refs.pathway_gene_link pgl
    INNER JOIN refs.ncbi_entrez n ON pgl.entrez_id = n.geneId
    {genes_string} GROUP BY path_id;""".format()

    size_dict = dict()
    try:
        con = mdb.connect(**dbvars)
        cur = con.cursor()
        cur.execute(cmd)
        # assert isinstance(cur.rowcount, int)
        row_count = cur.rowcount
        if not row_count:
            print "No pathways found."
            return size_dict
        for row_no in xrange(row_count):
            row = cur.fetchone()
            size_dict[row[0]] = row[1]

    except mdb.Error as e:
        print "Error %d: %s" % (e.args[0], e.args[1])
    finally:
        if con:
            con.close()
    return size_dict


def lookup_patient_counts(table_name, ignore_genes):
    """Get pathway sizes when exclude genes are provided."""

    if(ignore_genes):
        genes_string = repr(tuple(ignore_genes)).replace(",)", ")")
        genes_string = "WHERE hugo_symbol NOT IN {}".format(genes_string)
    else:
        genes_string = ""
    cmd = """SELECT patient_id, count(DISTINCT entrez_id)
              FROM {table_name} {genes_string} GROUP BY patient_id"""\
        .format(table_name=table_name, genes_string=genes_string)

    patient_size_dict = dict()
    try:
        con = mdb.connect(**dbvars)
        cur = con.cursor()
        cur.execute(cmd)
        # assert isinstance(cur.rowcount, int)
        row_count = cur.rowcount
        if not row_count:
            print "No pathways found."
            return patient_size_dict
        for row_no in xrange(row_count):
            row = cur.fetchone()
            patient_size_dict[row[0]] = row[1]

    except mdb.Error as e:
        print "Error %d: %s" % (e.args[0], e.args[1])
    finally:
        if con:
            con.close()
    return patient_size_dict


def build_path_patient_dict(table_name, ignore_genes):
    """Returns dict. Maps path_id (int) -> {Set of patient_ids}."""
    if ignore_genes:
        genes_string = repr(tuple(ignore_genes)).replace(",)", ")")
        genes_string = "WHERE hugo_symbol NOT IN {}".format(genes_string)
    else:
        genes_string = ""
    cmd = """SELECT pgl.path_id, patient_id FROM
            (SELECT DISTINCT patient_id, entrez_id FROM {table_name}
             {genes_string}) pg
            INNER JOIN
            refs.pathway_gene_link pgl ON pg.entrez_id = pgl.entrez_id;"""\
        .format(table_name=table_name, genes_string=genes_string)
    path_patient_dict = dict()
    try:
        con = mdb.connect(**dbvars)
        cur = con.cursor()
        cur.execute(cmd)
        # assert isinstance(cur.rowcount, int)
        row_count = cur.rowcount
        if not row_count:
            print "No patient-pathway pairs found."
            return path_patient_dict
        for row_no in xrange(row_count):
            row = cur.fetchone()
            path_id = row[0]
            patient_id = row[1]
            if path_id in path_patient_dict:
                path_patient_dict[path_id].add(patient_id)
            else:
                path_patient_dict[path_id] = {patient_id}
    except mdb.Error as e:
        print "Error %d: %s" % (e.args[0], e.args[1])
    finally:
        if con:
            con.close()
    return path_patient_dict


def fetch_path_ids_interest_genes(interest_genes):
    """Get pathway ids containing genes in (possibly empty) interest set."""
    rows = None
    all_path_ids = list()
    if interest_genes:
        genes_string = repr(tuple(interest_genes))
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


def get_gene_combs_hit(table_name):
        """Gets patient-pathway gene overlap info from databse.
        Only called by _populate_exclusive_cooccurring.
        """
        rows = None
        gene_lists = list()
        path_genes_dict = dict()

        cmd_maxlen = "SET group_concat_max_len = 10000;"
        cmd = """SELECT DISTINCT path_id, symbols FROM
            (
            # PATH, HUGO PAIRS in pathway of interest.
            SELECT path_id, group_concat(DISTINCT hugo_symbol
            ORDER BY hugo_symbol SEPARATOR ',') AS symbols
            FROM {table} t
            INNER JOIN refs.`pathway_gene_link` pgl
            ON t.entrez_id = pgl.entrez_id
            GROUP BY path_id, patient_id
            ) g;""". \
            format(table=table_name)
        try:
            con = mdb.connect(**dbvars)
            cur = con.cursor()
            cur.execute(cmd_maxlen)
            cur.execute(cmd)
            # rows = cur.fetchall()
            numrows = cur.rowcount
            for i in xrange(numrows):
                row = cur.fetchone()
                path_id = row[0]
                gene_list = row[1].split(',')
                if path_id in path_genes_dict:
                    path_genes_dict[path_id].append(gene_list)
                else:
                    path_genes_dict[path_id] = [gene_list]
        except mdb.Error as e:
            print "Error %d: %s" % (e.args[0], e.args[1])
        finally:
            if con:
                con.close()
        # prev returned list of lists, now dictionary with list of lists as vals
        return path_genes_dict


def get_gene_counts(table_name):
        """ Fetch dictionary of dictionaries: path -> gene -> patients with mutation.
        Dictionary may be empty if no pathway genes were mutated."""
        rows = None
        path_gene_dict = defaultdict(dict)
        cmd0 = """SET SESSION group_concat_max_len = 30000;"""
        # HUGO LIST AND PATIENT COUNTS
        cmd2 = """SELECT path_id, hugo_symbol, count(DISTINCT patient_id)
            AS n_patients, GROUP_CONCAT(DISTINCT patient_id) AS patients
            FROM {table} t
            # gene subset in pathway of interest
            INNER JOIN refs.`pathway_gene_link` pgl
            ON t.entrez_id = pgl.entrez_id
            GROUP BY path_id, hugo_symbol;""" \
            .format(table=table_name)
        try:
            con = mdb.connect(**dbvars)
            cur = con.cursor()
            cur.execute(cmd0)
            cur.execute(cmd2)
            row_count = cur.rowcount
            if not row_count:
                # NO GENES MUTATED. n_effective < n_pathway
                return path_gene_dict
            for i in xrange(row_count):
                row = cur.fetchone()
                path_id = row[0]
                gene = row[1]
                coverage = int(row[2])
                patient_names = row[3].split(',')
                if not len(patient_names) == coverage:
                    raise Exception(
                        "Pathway coverage query gives inconsistent " +
                        "patient counts and patient names; truncated "
                        "group_concat?")
                path_gene_dict[path_id][gene] = patient_names
        except mdb.Error as e:
            print "Error %d: %s" % (e.args[0], e.args[1])
        finally:
            if con:
                con.close()
        # OLD: count_dict : gene -> n_patients; total_patients
        # self.geneMatrix.add_gene_patients(gene, patient_names)
        return path_gene_dict