from collections import defaultdict, Counter
from . import db


def lookup_background_size(ignore_genes=None, alg=None, bmr_table=None):
    """Get background genome size sizes. Optional ignore_genes iterable.

    - Uses entrez_length table built from genes in msigdb (pathway_genes).
    - length_bp is cds length. or rna length for noncoding rna.
    - effective_bp is length_bp scaled by noncoding mutation rate (mutations/Mb)
        from mutsigcv paper, table s4 and s5. floor=1.5, ceiling=15.

    Args:
        ignore_genes (iterable): contains gene symbols to ignore.
        alg (str): specifies algorithm. ['gene_count', 'gene_length',
            or 'bmr_length']
        bmr_table: specify table name for custom bmr, else use
            entrez_length table. only applies to 'bmr_length' algorithm.

    Return:
        int: length. Number of bases for length, else number of genes.
    """
    if alg is None:
        raise Exception("Must specify algorithm type.")
    if alg not in ['gene_count', 'gene_length', 'bmr_length']:
        raise ValueError("Algorithm must be gene_count, gene_length, "
                         "or bmr_length.")
    if ignore_genes:
        genes_str = repr(tuple(ignore_genes)).replace(",)", ")")
        genes_str = "WHERE hugo_symbol NOT IN {}".format(genes_str)
    else:
        genes_str = ""
    table_name = bmr_table if bmr_table else 'refs.entrez_length'
    if alg == 'bmr_length':
        field_str = "sum(effective_bp)"
    elif alg == 'gene_length':
        field_str = "sum(length_bp)"
    else:
        field_str = "count(*)"
    cmd = """SELECT {field_str} FROM {table} {genes_str};""".format(
        field_str=field_str, table=table_name, genes_str=genes_str)

    result = db.session.execute(cmd)
    row_count = result.rowcount
    if row_count != 1:
        print("Background size lookup failed.")
        return
    row = result.fetchone()
    background_size = int(row[0])

    return background_size


def lookup_path_sizes(ignore_genes=None):
    """Get pathway sizes. Optional ignore_genes iterable."""
    if ignore_genes:
        genes_string = repr(tuple(ignore_genes)).replace(",)", ")")
        genes_string = "WHERE symbol NOT IN {}".format(genes_string)
    else:
        genes_string = ""
    cmd = """SELECT path_id, count(DISTINCT entrez_id)
    FROM refs.pathway_gene_link pgl
    INNER JOIN refs.ncbi_entrez n ON pgl.entrez_id = n.geneId
    {genes_string} GROUP BY path_id;""".format(genes_string=genes_string)
    size_dict = dict()

    result = db.session.execute(cmd)
    row_count = result.rowcount
    if not row_count:
        print("No pathways found.")
        return size_dict
    for row_no in range(row_count):
        row = result.fetchone()
        size_dict[row[0]] = row[1]
    return size_dict


def lookup_path_lengths(ignore_genes=None, alg=None, bmr_table=None):
    """Get bp length for each pathway, with optional exclude genes."""
    if alg is None:
        raise Exception("Must specify algorithm type.")
    if alg not in ['gene_length', 'bmr_length']:
        raise ValueError("Algorithm must be gene_length, or bmr_length.")
    if ignore_genes:
        genes_string = repr(tuple(ignore_genes)).replace(",)", ")")
        genes_string = "WHERE hugo_symbol NOT IN {}".format(genes_string)
    else:
        genes_string = ""
    table_name = bmr_table if bmr_table else 'refs.entrez_length'
    field_str = 'effective_bp' if alg == 'bmr_length' else 'length_bp'
    cmd = """SELECT path_id, sum({field_str}) AS bp FROM {table} l
        INNER JOIN refs.pathway_gene_link pgl ON l.entrez_id = pgl.entrez_id
        {genes_string}
        GROUP BY pgl.path_id;""".format(genes_string=genes_string,
                                        field_str=field_str, table=table_name)
    len_dict = dict()

    result = db.session.execute(cmd)
    row_count = result.rowcount
    if not row_count:
        print("No pathways found.")
        return len_dict
    for i in range(row_count):
        row = result.fetchone()
        len_dict[row[0]] = row[1]
    return len_dict


def lookup_patient_counts(table_name, ignore_genes):
    """Get patient gene counts."""

    if(ignore_genes):
        genes_string = repr(tuple(ignore_genes)).replace(",)", ")")
        genes_string = "WHERE hugo_symbol NOT IN {}".format(genes_string)
    else:
        genes_string = ""
    cmd = """SELECT patient_id, count(DISTINCT entrez_id)
              FROM {table_name} {genes_string} GROUP BY patient_id;"""\
        .format(table_name=table_name, genes_string=genes_string)

    patient_size_dict = dict()
    result = db.session.execute(cmd)
    row_count = result.rowcount
    if not row_count:
        print("No pathways found.")
        return patient_size_dict
    for row_no in range(row_count):
        row = result.fetchone()
        patient_size_dict[row[0]] = row[1]

    return patient_size_dict


def lookup_hypermutated_patients(table_name, cutoff=500):
    """Get patient ids for patients with >500 mutations (or specified cutoff)."""
    cmd = """SELECT patient_id FROM {table_name}
              GROUP BY patient_id HAVING count(*)>{cutoff};"""\
        .format(table_name=table_name, cutoff=cutoff)
    patient_list = []
    result = db.session.execute(cmd)
    for row in result:
        patient_list.append(row[0])
    return patient_list


def lookup_patient_lengths(table_name, ignore_genes):
    """Get patient bp lengths."""
    if(ignore_genes):
        genes_string = repr(tuple(ignore_genes)).replace(",)", ")")
        genes_string = "WHERE m.hugo_symbol NOT IN {}".format(genes_string)
    else:
        genes_string = ""
    cmd = """SELECT patient_id, count(*)
      FROM {table_name} m
      {genes_string} GROUP BY patient_id;"""\
        .format(table_name=table_name, genes_string=genes_string)
    patient_len_dict = dict()
    result = db.session.execute(cmd)
    for row in result:
        patient_len_dict[row[0]] = row[1]
    return patient_len_dict


def count_patients(table_name):
    """Get patient count for project."""
    cmd = """SELECT count(distinct patient_id) FROM {table_name};"""\
        .format(table_name=table_name)
    patient_count = None
    result = db.session.execute(cmd)
    row_count = result.rowcount
    if not row_count == 1:
        print("Non single result from patient count query.")
        return patient_count
    row = result.fetchone()
    patient_count = row[0]
    return patient_count


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
    result = db.session.execute(cmd)
    row_count = result.rowcount
    if not row_count:
        print("No patient-pathway pairs found.")
        return path_patient_dict
    for row in result:
        path_id = row[0]
        patient_id = row[1]
        if path_id in path_patient_dict:
            path_patient_dict[path_id].add(patient_id)
        else:
            path_patient_dict[path_id] = {patient_id}
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
    result = db.session.execute(cmd1)
    rowCount = result.rowcount
    if not rowCount:
        raise Exception(
            "Result contains %g rows Ids for pathway lookup." % rowCount)

    # rows is [[id,name],[id,name],...]
    for row in result:
        all_path_ids.append(int(row[0]))
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
    result = db.session.execute(cmd1)
    row_count = result.rowcount
    if not row_count > 1:
        raise Exception(
            "Result contains %g rows Ids for pathway lookup."
            % row_count)
    # rows is [[id,name],[id,name],...]
    for pair in result:
        path_id = int(pair[0])
        path_name = pair[1]
        pathway_dict[path_id] = path_name
    return pathway_dict

def get_pway_lenstats_dict(mutation_table, ignore_genes):
    """Get length stats for all mutated pathways."""
    rows = None
    pathway_lengths = dict()
    genes_string = ""
    if ignore_genes:
        genes_string = repr(tuple(ignore_genes)).replace(",)", ")")
        genes_string = "WHERE m.hugo_symbol NOT IN {}".format(genes_string)

    cmd1 = """SELECT g.path_id,
        cast(lmin/1000 AS DECIMAL(10,1)) AS `min`,
            group_concat(DISTINCT CASE e.`length_bp` WHEN lmin THEN m.hugo_symbol
                ELSE NULL END ORDER BY m.hugo_symbol SEPARATOR ', ') AS min_gene,
        cast(lmax/1000 AS DECIMAL(10,1)) AS `max`,
            group_concat(DISTINCT CASE e.`length_bp` WHEN lmax THEN m.hugo_symbol
                ELSE NULL END ORDER BY m.hugo_symbol SEPARATOR ', ') AS max_gene,
                lavg, lvar
         FROM (SELECT DISTINCT entrez_id, hugo_symbol FROM `{mutation_table}`) m
         INNER JOIN refs.entrez_length e ON m.entrez_id = e.entrez_id
        INNER JOIN refs.`pathway_gene_link` pgl ON e.entrez_id = pgl.entrez_id
        INNER JOIN # pway_stats
        (SELECT path_id,
        min(length_bp) AS `lmin`,
        max(length_bp) AS `lmax`,
        cast(AVG(length_bp)/1000 AS DECIMAL(10,1)) AS `lavg`,
        cast(var_samp(length_bp/1000) AS DECIMAL(10,1)) AS `lvar`
         FROM (SELECT DISTINCT hugo_symbol, entrez_id FROM `{mutation_table}`) m
         INNER JOIN refs.entrez_length e ON m.entrez_id = e.entrez_id
            INNER JOIN refs.`pathway_gene_link` pgl ON e.entrez_id = pgl.entrez_id
            {exclude_str} GROUP BY path_id) g ON g.path_id = pgl.`path_id`
            {exclude_str} GROUP BY g.path_id;"""\
        .format(mutation_table=mutation_table, exclude_str=genes_string)
    result = db.session.execute(cmd1)
    row_count = result.rowcount
    if not row_count > 1:
        raise Exception(
            "Result contains %g rows Ids for pathway lookup."
            % row_count)
    # rows is [[id,min,max,avg],[id,min,max,avg],...]
    for temp_lengths in result:
        path_id = int(temp_lengths[0])
        len_min = str(temp_lengths[1])
        gene_min = str(temp_lengths[2])
        len_max = str(temp_lengths[3])
        gene_max = str(temp_lengths[4])
        len_avg = str(temp_lengths[5])
        len_var = str(temp_lengths[6])
        pathway_lengths[path_id] = (len_min, gene_min, len_max, gene_max,
                                    len_avg, len_var)
    return pathway_lengths


def fetch_path_info_global():
    """Get url, brief description and contributor as tuple."""
    url_row = None
    cmd = "SELECT path_id, info_url, `description_brief`, contributor " \
          "FROM refs.pathways;"
    info_dict = dict()
    result = db.session.execute(cmd)
    row_count = result.rowcount
    if not row_count > 1:
        raise Exception("Failed info lookup for all pathways.")
    # rows is [[id,url,desc,contrib],[id,url,desc,contrib],...]
    for row in result:
        path_id = row[0]
        url = row[1]
        desc = row[2]
        contrib = row[3]
        info_dict[path_id] = dict(url=url, desc=desc, contrib=contrib)
    return info_dict


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
        db.session.execute(cmd_maxlen)
        result = db.session.execute(cmd)
        for row in result:
            path_id = row[0]
            gene_list = row[1].split(',')
            if path_id in path_genes_dict:
                path_genes_dict[path_id].append(gene_list)
            else:
                path_genes_dict[path_id] = [gene_list]
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
        db.session.execute(cmd0)
        result = db.session.execute(cmd2)
        row_count = result.rowcount
        if not row_count:
            # NO GENES MUTATED. n_effective < n_pathway
            return path_gene_dict
        for row in result:
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

        # OLD: count_dict : gene -> n_patients; total_patients
        # self.geneMatrix.add_gene_patients(gene, patient_names)
        return path_gene_dict


def get_annotation_dict(table_name):
    """ Fetch annotation dictionary: (hugo, patient) -> annot."""
    annot_dict = dict()
    cmd0 = """SET SESSION group_concat_max_len = 30000;"""
    # HUGO LIST AND PATIENT COUNTS
    cmd1 = """SELECT hugo_symbol, patient_id, GROUP_CONCAT(DISTINCT annot) AS annot
        FROM {table} t
        GROUP BY patient_id, hugo_symbol;""" \
        .format(table=table_name)

    db.session.execute(cmd0)
    result = db.session.execute(cmd1)
    for row in result:
        hugo, patient, annot = row
        annot_list = annot.split(',')
        # shorten annotation if necessary
        if len(annot_list) > 1:
            c = Counter(annot_list)  # Counter({'a': 5, 's': 1})
            annot = c.most_common()[0][0] + '+'
        annot_dict[(hugo, patient)] = annot
    return annot_dict
