import unicodedata as ud
from glob import glob
import os
import re
import random
import string
from itertools import combinations_with_replacement
import numpy as np
import subprocess
from collections import OrderedDict

import pandas as pd

class GeneListTester(object):
    """Holds regex string for testing gene lists."""
    hugo_re_str = "([A-Z0-9-]{1,}[,]{0,})+$"

    @staticmethod
    def is_valid(genes_str):
        if re.match(GeneListTester.hugo_re_str, genes_str):
            return True
        else:
            return False


def simplify_string(odd_string):
    keepcharacters = (' ', '.', '_', '-')
    safe_string = "".join(c for c in odd_string if c.isalnum()
                          or c in keepcharacters)\
        .rstrip().replace(' ', '_')
    if safe_string:
        safe_string = safe_string.encode('latin-1', 'ignore')\
            .decode('latin-1')
        safe_string = ud.normalize("NFC", safe_string)
    return safe_string


def zip_svgs(proj_dir, svgo_path='svgo'):
    """Zip and delete all svgs in matrix_svg, pathways_svg."""
    svg_dirs = ['pathways_svg', 'matrix_svg', '.', 'matrix_svg_cnv']
    for plot_dir in svg_dirs:
        plot_path = os.path.join(proj_dir, plot_dir)
        if not os.path.exists(plot_path):
            continue
        with open(os.devnull, "r") as fnullin:
            with open(os.devnull, "w") as fnullout:
                subprocess.check_call([svgo_path, '-f', plot_path],
                                      stdin=fnullin, stdout=fnullout,
                                      stderr=fnullout)
        file_names = glob(os.path.join(plot_path, '*.svg'))
        for file_name in file_names:
            # check_call(['/usr/local/bin/svgo', file_name])
            subprocess.check_call(['gzip', file_name])
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


def get_at_least_n(vals, n, nmax=200):
    """Get at least n values from list. Extend if following values are the same
    as nth value."""
    try:
        v_n = vals[n-1]
        v_next = vals[n]
    except IndexError:
        return vals[:n]

    if v_next == v_n and n < nmax:
        return get_at_least_n(vals, n+1, nmax)
    else:
        return vals[:n]

def strip_contributors(string_in):
    contribs = ['KEGG_', 'BIOCARTA_', 'ST_', 'SA_', 'SIG_', 'WNT_', 'PID_',
                'REACTOME_']
    for contrib in contribs:
        if string_in.startswith(contrib):
            return string_in[len(contrib):]
    return string_in


def get_nice_name(name_str):
    name_str = strip_contributors(name_str)
    name_str = name_str.replace('_', ' ')
    name_str = name_str.replace('TEL PATHWAY', 'TELOMERASE PATHWAY')
    name_str = name_str.replace('RNA PATHWAY', 'PKR SIGNALING PATHWAY')
    return name_str


def get_wait_time_string(wait_time):
    days = wait_time.days
    hours = wait_time.seconds/60/60
    mins = wait_time.seconds/60 - hours*60
    days_str = str(days) + 'day' if days else ''
    hours_str = str(hours) + 'h' if hours else ''
    mins_str = str(mins) + 'min' if mins else ''
    if days > 1:
        return str(days) + ' days'
    elif days:
        return ' '.join([days_str, hours_str]).strip()
    # < 1 day
    elif hours:
        return ' '.join([hours_str, mins_str]).strip()
    # < 1 hour
    elif mins:
        return mins_str
    # < 1 min
    else:
        return str(wait_time.seconds) + 's'


def generate_random_str(length=6):
    """ascii letters (upper+lower), numbers. http://stackoverflow.com/a/2257449"""
    charset = string.ascii_uppercase + string.digits + string.ascii_lowercase
    char_list = [random.SystemRandom().choice(charset) for _ in xrange(length)]
    return ''.join(char_list)


def objects_to_dataframe(object_list, attr_list, index_col=None,
                         header_map=dict(), dtype=None):
    """Create pandas dataframe from list of model instances.

    Args:
         attr_list (str[]): object attributes to turn into columns
         index_col (str): optional attribute to use as index
         header_map (dict): attribute name -> column name map
         dtype (type|dict) = optional dtype for dataframe. e.g. object
    """
    d = OrderedDict()
    for header in attr_list:
        new_header = header_map[header] if header in header_map else header
        d[new_header] = [getattr(obj, header) for obj in object_list]
    inds = [getattr(ob, index_col) for ob in object_list] if index_col else None
    df = pd.DataFrame(d, dtype=dtype, index=inds)
    df = df.where((pd.notnull(df)), None)  # None instead of nan
    return df
