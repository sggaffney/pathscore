from collections import defaultdict, Counter

from sqlalchemy import text, TextClause

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
        raise ValueError("Must specify algorithm type.")
    if alg not in ['gene_count', 'gene_length', 'bmr_length']:
        raise ValueError("Algorithm must be gene_count, gene_length, "
                         "or bmr_length.")
    exclude_genes_str = _get_gene_exclusion_sql(ignore_genes, symbol_col="hugo_symbol")
    table_name = bmr_table if bmr_table else 'refs.entrez_length'
    if alg == 'bmr_length':
        field_str = "sum(effective_bp)"
    elif alg == 'gene_length':
        field_str = "sum(length_bp)"
    else:
        field_str = "count(*)"
    cmd = text(f"""SELECT {field_str} FROM {table_name} {exclude_genes_str};""")
    result = db.session.execute(cmd)
    row_count = result.rowcount
    if row_count != 1:
        print("Background size lookup failed.")
        return
    row = result.fetchone()
    background_size = int(row[0])
    return background_size


def execute_cmd_fetchone_col0_col1_map(cmd: TextClause) -> dict:
    out_dict = dict()
    result = db.session.execute(cmd)
    row_count = result.rowcount
    if not row_count:
        print("No pathways found.")
        return out_dict
    for row_no in range(row_count):
        row = result.fetchone()
        out_dict[row[0]] = row[1]
    return out_dict


def lookup_path_sizes(ignore_genes=None) -> dict:
    """Get pathway sizes. Optional ignore_genes iterable."""
    exclude_genes_str = _get_gene_exclusion_sql(ignore_genes, symbol_col="symbol")
    cmd = text(f"""SELECT path_id, count(DISTINCT entrez_id)
    FROM refs.pathway_gene_link pgl
    INNER JOIN refs.ncbi_entrez n ON pgl.entrez_id = n.geneId
    {exclude_genes_str} GROUP BY path_id;""")
    return execute_cmd_fetchone_col0_col1_map(cmd)


def _get_gene_exclusion_sql(ignore_genes: list, symbol_col="hugo_symbol") -> str:
    if ignore_genes:
        genes_string = repr(tuple(ignore_genes)).replace(",)", ")")
        return f"WHERE {symbol_col} NOT IN {genes_string}"
    return ""


def _get_gene_inclusion_sql(require_genes: list, symbol_col="hugo_symbol") -> str:
    if require_genes:
        genes_string = repr(tuple(require_genes)).replace(",)", ")")
        return f"WHERE {symbol_col} IN {genes_string}"
    return ""


def lookup_path_lengths(ignore_genes=None, alg=None, bmr_table=None) -> dict:
    """Get bp length for each pathway, with optional exclude genes."""
    if alg is None:
        raise ValueError("Must specify algorithm type.")
    if alg not in ['gene_length', 'bmr_length']:
        raise ValueError("Algorithm must be gene_length, or bmr_length.")
    exclude_genes_str: str = _get_gene_exclusion_sql(ignore_genes)
    table_name = bmr_table if bmr_table else 'refs.entrez_length'
    field_str = 'effective_bp' if alg == 'bmr_length' else 'length_bp'
    cmd = text(f"""SELECT path_id, sum({field_str}) AS bp FROM {table_name} l
        INNER JOIN refs.pathway_gene_link pgl ON l.entrez_id = pgl.entrez_id
        {exclude_genes_str}
        GROUP BY pgl.path_id;""")
    return execute_cmd_fetchone_col0_col1_map(cmd)


def lookup_patient_counts(table_name: str, ignore_genes: list):
    """Get patient gene counts."""
    exclude_genes_string: str = _get_gene_exclusion_sql(ignore_genes)
    cmd = text(f"""SELECT patient_id, count(DISTINCT entrez_id)
               FROM {table_name} {exclude_genes_string} GROUP BY patient_id;""")
    return execute_cmd_fetchone_col0_col1_map(cmd)


def lookup_hypermutated_patients(table_name, cutoff=500):
    """Get patient ids for patients with >500 mutations (or specified cutoff)."""
    cmd = text(f"""SELECT patient_id FROM {table_name}
               GROUP BY patient_id HAVING count(*)>{cutoff};""")
    patient_list = []
    result = db.session.execute(cmd)
    for row in result:
        patient_list.append(row[0])
    return patient_list


def lookup_patient_lengths(table_name, ignore_genes):
    """Get patient bp lengths."""
    exclude_genes_str = _get_gene_exclusion_sql(ignore_genes, symbol_col="m.hugo_symbol")
    cmd = text(f"""SELECT patient_id, count(*)
        FROM {table_name} m
        {exclude_genes_str} GROUP BY patient_id;""")
    patient_len_dict = dict()
    result = db.session.execute(cmd)
    for row in result:
        patient_len_dict[row[0]] = row[1]
    return patient_len_dict


def count_patients(table_name):
    """Get patient count for project."""
    cmd = text(f"""SELECT count(distinct patient_id) FROM {table_name};""")
    patient_count = None
    result = db.session.execute(cmd)
    row_count = result.rowcount
    if not row_count == 1:
        print("Non single result from patient count query.")
        return patient_count
    row = result.fetchone()
    patient_count = row[0]
    return patient_count


def build_path_patient_dict(table_name, ignore_genes: list):
    """Returns dict. Maps path_id (int) -> {Set of patient_ids}."""
    exclude_genes_str = _get_gene_exclusion_sql(ignore_genes, symbol_col="hugo_symbol")
    cmd = text(f"""SELECT pgl.path_id, patient_id FROM
            (SELECT DISTINCT patient_id, entrez_id FROM {table_name}
             {exclude_genes_str}) pg
            INNER JOIN
            refs.pathway_gene_link pgl ON pg.entrez_id = pgl.entrez_id;""")
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
    all_path_ids = list()
    require_genes_str = _get_gene_inclusion_sql(interest_genes, symbol_col="symbol")
    cmd1 = text(f"""SELECT distinct path_id FROM refs.pathway_gene_link pgl
        INNER JOIN (SELECT geneId FROM refs.ncbi_entrez
            {require_genes_str}) g
        ON pgl.entrez_id = g.geneId ORDER BY path_id;""")
    result = db.session.execute(cmd1)
    row_count = result.rowcount
    if not row_count:
        raise Exception(
            "Result contains %g rows Ids for pathway lookup." % row_count)
    # result is [[id,name],[id,name],...]
    for row in result:
        all_path_ids.append(int(row[0]))
    return all_path_ids


def get_pathway_name_dict():
    """Gets name for all pathways, stored in dict: pathid -> pathname."""
    pathway_dict = dict()
    # GET pway_size
    cmd1 = text("""SELECT p.path_id, pathway_name FROM refs.pathways p
        INNER JOIN
        (SELECT DISTINCT path_id FROM refs.pathway_gene_link) l
        ON p.path_id = l.path_id;""")
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
    pathway_lengths = dict()
    exclude_genes_str = _get_gene_exclusion_sql(ignore_genes, symbol_col="m.hugo_symbol")
    cmd1 = text(f"""SELECT g.path_id,
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
            {exclude_genes_str} GROUP BY path_id) g ON g.path_id = pgl.`path_id`
            {exclude_genes_str} GROUP BY g.path_id;""")
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
    cmd = text("SELECT path_id, info_url, `description_brief`, contributor "
               "FROM refs.pathways;")
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
    cmd_maxlen = text("SET group_concat_max_len = 10000;")
    cmd = text(f"""SELECT DISTINCT path_id, symbols FROM
        (
        # PATH, HUGO PAIRS in pathway of interest.
        SELECT path_id, group_concat(DISTINCT hugo_symbol
        ORDER BY hugo_symbol SEPARATOR ',') AS symbols
        FROM {table_name} t
        INNER JOIN refs.`pathway_gene_link` pgl
        ON t.entrez_id = pgl.entrez_id
        GROUP BY path_id, patient_id
        ) g;""")
    db.session.execute(cmd_maxlen)
    path_genes_dict = dict()
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
    cmd0 = text("""SET SESSION group_concat_max_len = 30000;""")
    # HUGO LIST AND PATIENT COUNTS
    cmd2 = text(f"""SELECT path_id, hugo_symbol, count(DISTINCT patient_id)
        AS n_patients, GROUP_CONCAT(DISTINCT patient_id) AS patients
        FROM {table_name} t
        # gene subset in pathway of interest
        INNER JOIN refs.`pathway_gene_link` pgl
        ON t.entrez_id = pgl.entrez_id
        GROUP BY path_id, hugo_symbol;""")
    db.session.execute(cmd0)
    result = db.session.execute(cmd2)
    row_count = result.rowcount
    path_gene_dict = defaultdict(dict)
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
    """Fetch annotation dictionary: (hugo, patient) -> annot."""
    cmd0 = text("""SET SESSION group_concat_max_len = 30000;""")
    # HUGO LIST AND PATIENT COUNTS
    cmd1 = text(f"""SELECT hugo_symbol, patient_id, GROUP_CONCAT(DISTINCT annot) AS annot
        FROM {table_name} t
        GROUP BY patient_id, hugo_symbol;""")
    annot_dict = dict()
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
