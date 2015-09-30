from flask import current_app
import os


def get_user_folder(user_id):
    data_root = current_app.config['DATA_ROOT']
    return os.path.join(data_root, str(user_id))


def get_project_folder(upload_obj):
    user_folder = get_user_folder(upload_obj.user_id)
    return os.path.join(user_folder, str(upload_obj.file_id))


def get_js_name(upload_obj):
    return upload_obj.get_local_filename() + ".js"


def get_js_path(upload_obj):
    dir_path = get_project_folder(upload_obj)
    js_name = get_js_name(upload_obj)
    return os.path.join(dir_path, js_name)


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

