import os
import errno
from itertools import product

import pandas as pd
import numpy as np

from scipy.spatial.distance import pdist, squareform
from sklearn import manifold
# from sklearn.metrics import jaccard_similarity_score, euclidean_distances, \
#     hamming_loss
from sklearn.decomposition import PCA

import app
from . import naming_rules, db
from .models import UserFile
from . import get_effective_pathways as gep
from . import db_lookups as lookup


metric_list = ['jaccard', 'dice', 'odds_inv', 'hamming', 'rogerstanimoto',
               'sokalsneath']  # yule
alg_list = ['MDS', 'NMDS']


def get_table_names():
    result = app.db.engine.execute('show tables;')
    table_names = [i[0] for i in result.fetchall()]
    result.close()
    return table_names


def load_proj_paths(user_upload):
    """Build dictionary of project attributes for MDS plot."""
    # scores_path, names_path, tree_path = naming_rules.get_tree_score_paths(
    #     user_upload)
    detail_path = naming_rules.get_detailed_path(user_upload)
    all_paths = gep.load_pathway_list_from_file(detail_path)
    # build table for loading lenstats
    proj_dir = naming_rules.get_project_folder(user_upload)
    data_path = os.path.join(proj_dir, user_upload.get_local_filename())
    table_name = user_upload.get_table_name()
    create_table = table_name not in get_table_names()

    if create_table:
        gep.MutationTable(table_name, data_path,
                          unused_path=None,
                          rejected_path=None,
                          has_annot=user_upload.has_annot)
    ignore_genes = str(user_upload.ignore_genes).split(
        ',') if user_upload.ignore_genes else []
    # alg = user_upload.algorithm
    # GENOME SIZE
    # genome_size = lookup.lookup_background_size(ignore_genes=ignore_genes,
    #                                             alg=alg, bmr_table=None)
    # HYPERMUTATED
    hypermutated = lookup.lookup_hypermutated_patients(table_name)
    # Pathway gene lists and patient-gene pairs
    #  e.g. (1035L, [[u'PLCB1']]); (1035L, {u'PLCB1': [u'71347']})
    #  path_genelists_dict = g.get_gene_combs_hit(table_name)
    path_patient_dict = lookup.build_path_patient_dict(table_name, ignore_genes)
    patient_length_dict = lookup.lookup_patient_lengths(table_name,
                                                        ignore_genes)
    # path_genepatients_dict = lookup.get_gene_counts(table_name)
    # n_patients = user_upload.n_patients

    patient_list = patient_length_dict.keys()
    proj_name = user_upload.get_local_filename()

    proj_info = dict(proj_name=proj_name,
                     all_paths=all_paths,
                     path_patient_dict=path_patient_dict,
                     patient_list=patient_list,
                     hypermutated=hypermutated,
                     proj_dir=proj_dir)
    if create_table:
        gep.drop_table(table_name)
    return proj_info


def build_mds_dataframe(proj_info):
    """Create dataframe to hold pathway-patient is-mutated boolean values.

    :param proj_info: dict
    :return: pd.DataFrame
    """
    hypermutated = proj_info.get('hypermutated', None)
    path_patient_dict = proj_info['path_patient_dict']
    patient_list = proj_info['patient_list']
    all_paths = proj_info['all_paths']
    all_path_ids = path_patient_dict.keys()
    df = pd.DataFrame(index=all_path_ids, columns=patient_list, dtype=bool,
                      data=False)
    for path_id, patients in path_patient_dict.items():
        df.loc[path_id, list(patients)] = True
    use_paths = [i for i in all_paths if i.gene_coverage]
    use_paths_ids = [int(i.path_id) for i in use_paths]
    non_hypermutated = [i for i in patient_list if i not in hypermutated]
    df = df.loc[use_paths_ids, sorted(non_hypermutated)].copy()
    return df


class EchoDict(dict):
    def __missing__(self, key):
        return key


def get_inv_odds_ratio(u, v):
    a = u.astype(np.bool)
    b = v.astype(np.bool)
    return np.float64((a & ~b).sum() + 1) * \
        ((~a & b).sum() + 1) / ((a & b).sum() + 1) / ((~a & ~b).sum() + 1)


metric_dict = EchoDict(odds_inv=get_inv_odds_ratio)


def mkdir_p(path):
    """Replicates mkdir -p functionality.

    via https://stackoverflow.com/a/600612"""
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def build_mds_pts(df, proj_dir='', make_plot=False):
    """Create numpy files holding mds points for multiple metrics/algorithms."""
    out_dir = os.path.join(proj_dir, 'npy_pts')
    if ~os.path.exists(out_dir):
        mkdir_p(out_dir)
    # iterate throguh each metric / alg pair
    for metric, mds_alg in product(metric_list, alg_list):
        x = np.array(df)
        sim = squareform(pdist(x, metric_dict[metric]))
        seed = np.random.RandomState(seed=2)
        if mds_alg == 'NMDS':
            col = 'b'
            mds = manifold.MDS(n_components=2, metric=False, max_iter=3000,
                               eps=1e-12,
                               dissimilarity="precomputed", random_state=seed,
                               n_jobs=1,
                               n_init=1)
        else:
            mds = manifold.MDS(n_components=3, max_iter=3000, eps=1e-12,
                               random_state=seed,
                               dissimilarity="precomputed", n_jobs=1)
            col = 'r'
        base_str = '{}_{}'.format(metric, mds_alg)
        pos = mds.fit(sim).embedding_
        # RESCALE: pos *= np.sqrt((X ** 2).sum()) / np.sqrt((pos ** 2).sum())
        # ROTATE
        clf = PCA(n_components=2)
        pos = clf.fit_transform(pos)
        out_path_npy = os.path.join(out_dir, '{}.npy'.format(base_str))
        np.save(out_path_npy, pos)

        if make_plot:
            import matplotlib.pyplot as plt
            hfig, ax = plt.subplots(1, figsize=(5, 5))
            plt.scatter(pos[:, 0], pos[:, 1], lw=0, label='MDS',
                        color=col)  # s=50, color='turquoise'
            plt.title(base_str)
            # plt.legend(scatterpoints=1, loc='best', shadow=False)
            out_path_png = os.path.join(out_dir,
                                        'scatter_{}.png'.format(base_str))
            hfig.savefig(out_path_png)


def build_proj_mds(proj_id):
    user_upload = UserFile.query.get(proj_id)  # type: UserFile
    proj_info = load_proj_paths(user_upload)
    df = build_mds_dataframe(proj_info)
    mkdir_p(naming_rules.get_mds_pts_dir(user_upload))
    tsv_path = naming_rules.get_mds_df_path(user_upload)
    df.to_csv(tsv_path, sep='\t')
    proj_dir = naming_rules.get_project_folder(user_upload)
    build_mds_pts(df, proj_dir=proj_dir, make_plot=False)
    user_upload.has_mds = True
    db.session.add(user_upload)
    db.session.commit()
    print("Successfully created MDS dataframe and metric-specific points.")
