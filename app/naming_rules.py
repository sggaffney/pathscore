from flask import current_app
import os


def get_user_folder(user_id):
    upload_folder = current_app.config['UPLOAD_FOLDER']
    return os.path.join(upload_folder, str(user_id))


def get_project_folder(upload_obj):
    user_folder = get_user_folder(upload_obj.user_id)
    return os.path.join(user_folder, str(upload_obj.file_id))


def get_tree_score_paths(upload_obj):
    """Get path tuple (scores_path, scored_names) from upload object."""
    proj_suffix = upload_obj.get_local_filename()
    dir_path = get_project_folder(upload_obj)
    scores_path = os.path.join(dir_path, proj_suffix + '_scores.txt')
    scored_names = os.path.join(dir_path, proj_suffix + '_score_names.txt')
    svg_out_path = os.path.join(dir_path, proj_suffix + '_tree.svg')
    return scores_path, scored_names, svg_out_path


def get_detailed_path(upload_obj):
    """Get detailed path for upload object."""
    dir_path = get_project_folder(upload_obj)
    base_str = 'pathways_pvalues_'
    proj_suffix = upload_obj.get_local_filename()
    name_postfix = '_detail.txt'
    return os.path.join(dir_path, base_str + proj_suffix + name_postfix)


def get_apache_path(full_path):
    return full_path.replace('/www/pway', '/static/data')

