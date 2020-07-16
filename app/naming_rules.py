from flask import current_app
import os


mds_dir = 'npy_pts'


def get_user_folder(user_id):
    data_root = current_app.config['DATA_ROOT']
    return os.path.join(data_root, str(user_id))


def get_project_folder(upload_obj):
    user_folder = get_user_folder(upload_obj.user_id)
    return os.path.join(user_folder, str(upload_obj.file_id))


def get_bmr_folder(user_id):
    user_folder = get_user_folder(user_id)
    return os.path.join(user_folder, 'bmr')


def get_user_comparison_dir(user_id):
    """Holds directories with specific project list comparisons."""
    dir_path = get_user_folder(user_id)
    comparison_dir = os.path.join(dir_path, 'comparisons')
    return comparison_dir


def get_comparison_dir(user_id, proj_id_list):
    """Directory for specific project list comparison tables."""
    parent = get_user_comparison_dir(user_id)
    proj_str = '_'.join([str(i) for i in proj_id_list])
    path = os.path.join(parent, proj_str)
    return path


def get_js_name(upload_obj):
    return upload_obj.get_local_filename() + ".js"


def get_js_path(upload_obj):
    dir_path = get_project_folder(upload_obj)
    js_name = get_js_name(upload_obj)
    return os.path.join(dir_path, js_name)


def get_params_path(upload_obj):
    dir_path = get_project_folder(upload_obj)
    return os.path.join(dir_path,
                        upload_obj.get_local_filename() + "_params.txt")


def get_pval_path(upload_obj):
    dir_path = get_project_folder(upload_obj)
    root_name = upload_obj.get_local_filename()
    pval_name = root_name + '_pvalues.txt'
    return os.path.join(dir_path, pval_name)


def get_all_pways_path(upload_obj):
    dir_path = get_project_folder(upload_obj)
    root_name = upload_obj.get_local_filename()
    base_name = root_name + '_all_pways.tsv'
    return os.path.join(dir_path, base_name)


def get_tree_score_paths(upload_obj):
    """Get path tuple (scores_path, scored_names) from upload object."""
    proj_suffix = upload_obj.get_local_filename()
    dir_path = get_project_folder(upload_obj)
    scores_path = os.path.join(dir_path, proj_suffix + '_scores.txt')
    scored_names = os.path.join(dir_path, proj_suffix + '_score_names.txt')
    svg_out_path = os.path.join(dir_path, proj_suffix + '_tree.svg')
    return scores_path, scored_names, svg_out_path


def _get_pvalue_root(upload_obj):
    """Get detailed path for upload object."""
    dir_path = get_project_folder(upload_obj)
    base_str = 'pathways_pvalues_'
    proj_suffix = upload_obj.get_local_filename()
    return os.path.join(dir_path, base_str + proj_suffix)


def get_pvalue_path(upload_obj):
    return _get_pvalue_root(upload_obj) + '.txt'


def get_detailed_path(upload_obj):
    """Get detailed path for upload object."""
    name_prefix = _get_pvalue_root(upload_obj)
    name_postfix = '_detail.txt'
    return name_prefix + name_postfix


def get_hypermutated_path(upload_obj):
    hyper_file = 'hypermutated_{}.txt'.format(upload_obj.get_local_filename())
    dir_path = get_project_folder(upload_obj)
    return os.path.join(dir_path, hyper_file)


def get_target_dir(upload_obj):
    """Target plot directory. SEE ALSO get_target_path."""
    dir_path = get_project_folder(upload_obj)
    svg_dir = os.path.join(dir_path, 'pathways_svg')
    return svg_dir


def get_target_path(upload_obj, path_id, compressed=False):
    """Target plot path. Set compressed to True for svgz extension."""
    svg_dir = get_target_dir(upload_obj)
    ext = 'svg' if not compressed else 'svgz'
    svg_path = os.path.join(svg_dir, '{}.{}'.format(path_id, ext))
    return svg_path


def get_matrix_dir(upload_obj):
    """Target plot directory. SEE ALSO get_target_path."""
    dir_path = get_project_folder(upload_obj)
    svg_dir = os.path.join(dir_path, 'matrix_svg')
    return svg_dir


def get_matrix_path(upload_obj, path_id, compressed=False):
    """Target plot path. Set compressed to True for svgz extension."""
    svg_dir = get_matrix_dir(upload_obj)
    ext = 'svg' if not compressed else 'svgz'
    svg_path = os.path.join(svg_dir, '{}.{}'.format(path_id, ext))
    return svg_path


def get_mds_df_path(upload_obj):
    """Target plot path. Set compressed to True for svgz extension."""
    dir_path = get_project_folder(upload_obj)
    proj_suffix = upload_obj.get_local_filename()
    out_name = "{}_mds.tsv".format(proj_suffix)
    df_path = os.path.join(dir_path, mds_dir, out_name)
    return df_path


def get_mds_pts_dir(upload_obj):
    """Target plot path. Set compressed to True for svgz extension."""
    dir_path = get_project_folder(upload_obj)
    pts_path = os.path.join(dir_path, 'npy_pts')
    return pts_path


def get_unused_gene_path(upload_obj):
    unused_file = '{}_unused.txt'.format(upload_obj.get_local_filename())
    dir_path = get_project_folder(upload_obj)
    return os.path.join(dir_path, unused_file)


def get_rejected_gene_path(upload_obj):
    rejected_file = '{}_rejected.txt'.format(upload_obj.get_local_filename())
    dir_path = get_project_folder(upload_obj)
    return os.path.join(dir_path, rejected_file)


def get_apache_path(full_path):
    return full_path.replace('/www/pway', '/static/data')
