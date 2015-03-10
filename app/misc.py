import unicodedata as ud
from subprocess import check_call
from glob import glob
import os
from itertools import combinations_with_replacement
import numpy as np


def simplify_string(odd_string):
    keepcharacters = (' ', '.', '_')
    safe_string = "".join(c for c in odd_string if c.isalnum()
                          or c in keepcharacters)\
        .rstrip().replace(' ', '_')
    if safe_string:
        safe_string = safe_string.encode('latin-1', 'ignore')\
            .decode('latin-1')
        safe_string = ud.normalize("NFC", safe_string)
    return safe_string


def zip_svgs(proj_dir):
    """Zip and delete all svgs in matrix_svg, pathways_svg."""
    svg_dirs = ['pathways_svg', 'matrix_svg', '.']
    for plot_dir in svg_dirs:
        file_names = glob(os.path.join("{}".format(proj_dir), plot_dir, '*.svg'))
        for file_name in file_names:
            # check_call(['/usr/local/bin/svgo', file_name])
            check_call(['gzip', file_name])
            os.rename(file_name + '.gz', file_name + 'z')


def html_quotes(str_before):
    """Replaces apostrophe with html for single right apostrophe (&#8217;)."""
    return str_before.replace("'","&#39;").replace('"',"&quot;")


def get_pway_dissimilarity(gene_set1, gene_set2):
    """Takes two gene sets and calculates dissimilarity score from overlap."""
    overlap = len(gene_set1.intersection(gene_set2))
    min_len = min(len(gene_set1), len(gene_set2))
    return 1-float(overlap)/min_len


def get_distance_matrix(gene_sets):
    """Return symmetric matrix of dissimilarities from list of gene sets."""
    n_pways = len(gene_sets)
    scores = np.zeros([n_pways, n_pways])
    for ind_pair in combinations_with_replacement(xrange(n_pways), 2):
        score = get_pway_dissimilarity(set(gene_sets[ind_pair[0]]),
                                       set(gene_sets[ind_pair[1]]))
        scores[ind_pair[0], ind_pair[1]] = score
    # turn diagonal matrix into symmetric matrix
    scores = scores + np.transpose(scores)
    return scores


def get_at_least_n(vals, n):
    """Get at least n values from list. Extend if following values are the same
    as nth value."""
    try:
        v_n = vals[n-1]
        v_next = vals[n]
    except IndexError:
        return vals[:n]

    if v_next == v_n:
        return get_at_least_n(vals, n+1)
    else:
        return vals[:n]
