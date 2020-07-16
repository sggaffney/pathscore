import re
import logging
from collections import OrderedDict, Counter

import pandas as pd
import scipy.stats as stats


_logger = logging.getLogger(__name__)


STRIP_PREFIXES = {'KEGG_', 'REACTOME_', 'PID_', 'BIOCARTA_'}
# GREENS = sns.light_palette("green", n_colors=150).as_hex()
GREENS = ['#e5ffe5', '#e5ffe5', '#e3fee3', '#e1fde1', '#e0fce0', '#defbde', '#dcfadc', '#dbf9db', '#daf9da', '#d8f8d8', '#d6f7d6', '#d5f6d5', '#d3f5d3', '#d2f4d2', '#d0f3d0', '#cff3cf', '#cdf2cd', '#cbf1cb', '#caf0ca', '#c9efc9', '#c7eec7', '#c5edc5', '#c4edc4', '#c2ecc2', '#c1ebc1', '#c0eac0', '#bee9be', '#bce8bc', '#bae7ba', '#b9e7b9', '#b8e6b8', '#b6e5b6', '#b5e4b5', '#b3e3b3', '#b1e2b1', '#b0e1b0', '#afe1af', '#ade0ad', '#abdfab', '#a9dea9', '#a8dda8', '#a6dca6', '#a5dba5', '#a4dba4', '#a2daa2', '#a0d9a0', '#9ed89e', '#9ed79e', '#9cd69c', '#9ad59a', '#99d599', '#97d497', '#95d395', '#94d294', '#93d193', '#91d091', '#8fcf8f', '#8ecf8e', '#8cce8c', '#8bcd8b', '#89cc89', '#88cb88', '#86ca86', '#84c984', '#83c983', '#82c882', '#80c780', '#7ec67e', '#7dc57d', '#7bc47b', '#7ac37a', '#79c379', '#77c277', '#75c175', '#73c073', '#72bf72', '#70be70', '#6fbd6f', '#6dbc6d', '#6cbc6c', '#6abb6a', '#68ba68', '#67b967', '#66b866', '#64b764', '#62b662', '#61b661', '#5fb55f', '#5eb45e', '#5db35d', '#5bb25b', '#59b159', '#57b057', '#56b056', '#55af55', '#53ae53', '#52ad52', '#50ac50', '#4eab4e', '#4daa4d', '#4caa4c', '#4aa94a', '#48a848', '#47a747', '#45a645', '#44a544', '#42a442', '#41a441', '#3fa33f', '#3da23d', '#3ca13c', '#3aa03a', '#399f39', '#379e37', '#369e36', '#349d34', '#329c32', '#319b31', '#309a30', '#2e992e', '#2c982c', '#2b982b', '#299729', '#289628', '#269526', '#259425', '#239323', '#219221', '#209220', '#1f911f', '#1d901d', '#1b8f1b', '#1a8e1a', '#188d18', '#178c17', '#168c16', '#148b14', '#128a12', '#108910', '#0f880f', '#0e870e', '#0c860c', '#0b860b', '#098509', '#078407', '#058305', '#058205', '#038103', '#018001', '#008000']


def parse_cov_struct(cov_struct):
    """Convert PathScore coverage string to coverage dictionary.

    Example input: "struct('PTPN13',20.00,'RGS3',20.00,'SRC',20.00)"
    """
    try:
        vals = cov_struct[7:][:-1].replace("'", "").split(',')
    except TypeError:
        return {}
    genes = vals[0::2]
    covs = [float(i) for i in vals[1::2]]
    return dict(zip(genes, covs))


def gather_pathways(categ_names, tsv_dict, n_dict):
    """Build table for PathScore project comparison.

    Args:
        categ_names (iterable): Project names
        tsv_dict (dict): maps Project name -> tsv path for all_pways file.
        n_dict (dict): Project name -> number of cases in project. Note: should
            also include cases without any variants.
    """
    df_dict = {}
    for categ in categ_names:
        n_samples = n_dict[categ]
        df = pd.read_csv(tsv_dict[categ], sep='\t')
        df.rename(columns={'n_geness_total': 'n_genes',
                           'n_genes_total': 'n_genes'}, inplace=True)
        df['n_miss'] = n_samples - df['n_cov']
        df['rank'] = df.index.values + 1
        df = df.set_index(['pathway_id', 'pathway_name'])
        df.columns = pd.MultiIndex.from_product([[categ], df.columns])
        df_dict[categ] = df.copy()

    df = pd.concat([df_dict[i] for i in categ_names], axis=1)
    df.columns = df.columns.reorder_levels([1, 0])

    df['pc_cov'] = df['pc_cov'].fillna(0).astype(int)
    df['n_cov'] = df['n_cov'].fillna(0).astype(int)
    df['n_genes_mutated'] = df['n_genes_mutated'].fillna(0).astype(int)
    df[('n_genes', '')] = df['n_genes'].apply(
        lambda v: v.dropna().iloc[0], axis=1).astype(int)
    df.drop(columns=[('n_genes', v) for v in categ_names], inplace=True)
    df['gene_coverage'] = df['gene_coverage'].applymap(parse_cov_struct)
    # Add condensed gene list
    genes = df['gene_coverage'].applymap(lambda v: _get_genes_display(v))
    genes.columns = pd.MultiIndex.from_tuples([('genes', i) for i in genes.columns])

    pcn = df['pc_cov'].applymap(lambda v: f'{v}%') + df['n_cov'].applymap(lambda v: f' ({v})')
    pcn.columns = pd.MultiIndex.from_tuples([('pc_n', i) for i in pcn.columns])

    df = pd.concat([df, genes, pcn], axis=1)
    # Sort columns
    df = df[df.columns.sortlevel(0)[0]].copy()

    for categ in categ_names:
        df[('n_miss', categ)] = df['n_miss'][categ].fillna(n_dict[categ]).astype(int)

    df['rank'] = df['rank'].fillna(df['rank'].max().max()).astype(int)
    # n_cov, n_miss, pc_cov, n_genes_mutated, n_genes, p_value, rank
    temp = df['n_effective'] / df['n_actual']
    for categ in categ_names:
        df[('effect_floor', categ)] = temp[categ].where(temp[categ].fillna(0) > 1, 1)
    # Get shortened pathway names without database
    short_names = get_short_names(df.index.levels[1].values)
    df.index.set_levels(short_names, level='pathway_name', inplace=True)
    return df


def identify_enrichment(ref_name, compared_name, df, other_compared=None,
                        excel_out=False,
                        ):
    """Compare enrichment for project pair.

    Args:
        ref_name (str): Project name for enriched pathways
        compared_name (str): Project name for comparison
        df (pd.DataFrame): pathway stats table, output from gather_pathways.
        other_compared (iterable): other project names for calculating
            proportion P value.
        excel_out (bool): if True, write comparisons to Excel files.
    Returns:
        dictionary of project name -> enrichment dataframe
    """

    show_cols = ['pc_n', 'p_fisher', 'n_genes_mutated', 'n_genes', 'mut_genes',
                 'p_value', 'effect_floor', 'rank', ]  # 'pc_cov', 'n_cov', n_miss

    all_names = [i for i in df.columns.levels[1] if i]
    names_rest = [i for i in all_names if i != ref_name]
    other_compared = [] if other_compared is None else other_compared

    _logger.info(f"Comparing enrichment in {ref_name} to {compared_name}.")
    temp = df.copy()
    temp['mut_genes'] = temp[('genes', ref_name)]
    temp['effect_floor_diff'] = \
        temp['effect_floor'][ref_name] - temp['effect_floor'][compared_name]
    temp['p_fisher'] = temp.apply(
        lambda r: stats.fisher_exact([r['n_cov'][[ref_name, compared_name]],
                                      r['n_miss'][[ref_name, compared_name]]],
                                     alternative='greater')[1], axis=1)
    for categ in other_compared:
        temp[('p_fisher', categ)] = temp.apply(
            lambda r: stats.fisher_exact([r['n_cov'][[ref_name, categ]],
                                          r['n_miss'][[ref_name, categ]]],
                                         alternative='greater')[1], axis=1)
        # temp[('p_fisher', 'min')] = temp['p_fisher'].min(axis=1)
    # temp.sort_values('effect_floor_diff', ascending=False, inplace=True)
    temp['low_p'] = temp['p_value'][ref_name] < 0.05
    # temp[low_p]

    temp.sort_values(('rank', ref_name), ascending=True, inplace=True)
    out = temp.loc[(temp['low_p'] & (temp['effect_floor_diff'] > 0))][show_cols]
    out = out.reset_index(level=1)
    out['pathway_name'] = out['pathway_name'].transform(
        lambda v: v.replace('_', ' '))
    if excel_out:
        out.to_excel(f'enrichment_{ref_name}_gt_{compared_name}.xlsx')
    return out.copy()


def get_styled_comparison(comp, caption):
    """Style enrichment comparison table.

    Args:
        comp (pd.DataFrame): comparison table via identify_enrichment.
        caption (str): caption text for output.

    Returns:
         Styler object.
    """
    # Make rows yellow for low p_fisher
    rename_dict = {'pathway_name': caption,
                   'pc_n': 'Cases with variant',
                   'n_genes': 'Number of genes in pathway',
                   'n_genes_mutated': 'Number of genes mutated in pathway',
                   'p_fisher': 'P value, case proportion comparison',
                   'p_value': 'P value for pathway enrichment',
                   'rank': 'Pathway rank, by effect size',
                   'mut_genes': 'Pathway genes with mutation',
                   'effect_floor': 'PathScore effect size (floor=1)',
                   }
    comp = comp.rename(columns=rename_dict)
    styler = comp.style.apply(_style_using_p_values, axis=1)
    # Color coverage by percentage
    styler = styler.applymap(_style_coverage, subset=['Cases with variant'])
    # number formats
    styler = styler.format('{:.2f}', subset=['P value for pathway enrichment',
                                             'P value, case proportion comparison']).\
        format('{:.1f}', subset=['PathScore effect size (floor=1)'])
    # effect_max = comp.effect_floor.max().max()
    # styler.bar(subset=['effect_floor'], color='#5fba7d', vmin=1, vmax=10)
    styler.set_caption(caption)
    return styler


def _style_coverage(pcn_val):
    """Extract percentage from string like '82% (9)' and apply colormap."""
    pc = int(re.match('(\d+)%', pcn_val).groups()[0])
    hex_color = GREENS[pc]
    return f"background-color: {hex_color}"


def _style_using_p_values(s, fisher_cutoff=0.05, p_cutoff=0.05):
    vals = ['' for i in range(len(s))]
    # yellow rows for low fisher p
    if s[('P value, case proportion comparison', '')] < fisher_cutoff:
        vals = ['background-color: yellow' for i in range(len(s))]
        # make first column (pathway name) bold
        vals[0] += '; font-weight: bold'
        # fisher_ind = list(s.index.get_loc('p_fisher')).index(True)
    # red color for p values below 0.05
    for ind, col in enumerate(s.index.get_level_values(0)):
        if col.startswith('P value'):
            vals[ind] += '; number-format: 0.00'
            if s.iloc[ind] < p_cutoff:
                vals[ind] += '; color: red; font-weight: bold'
        if 'effect' in col.lower() and 'rank' not in col.lower():
            if s.iloc[ind] != 1:
                vals[ind] += '; number-format: 0.0'
    return vals


def _get_genes_display(cov_dict, n_show=4):
    keys = list(cov_dict.keys())
    genes = keys[:n_show]
    if len(keys) > n_show:
        genes += ['...']
    genes_show = ' '.join(genes)
    return genes_show


def get_short_names(pathway_names):
    """Get pathway names without database prefixes, excluding duplicates."""
    short_names = [_get_sm_name(i) for i in pathway_names]
    rep_names = set()
    for short_name, n in Counter(short_names).most_common():
        if n > 1:
            rep_names.add(short_name)
        else:
            break
    short_names = [_get_sm_name_reps(i, rep_names) for i in pathway_names]
    return short_names


def _get_sm_name(long_name):
    """Used by get_short_names."""
    for i in STRIP_PREFIXES:
        if long_name.startswith(i):
            return long_name[len(i):]
    return long_name


def _get_sm_name_reps(long_name, rep_names):
    """Used by get_short_names."""
    for i in STRIP_PREFIXES:
        if long_name.startswith(i):
            db = i[:-1]
            short_name = long_name[len(i):]
            if short_name in rep_names:
                short_name = '_'.join([short_name, db])
            return short_name
    return long_name
