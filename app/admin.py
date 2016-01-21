from datetime import datetime, timedelta
from models import UserFile, User, Role
from flask import current_app, render_template
from . import db
import os
import shutil
from threading import Thread
import time
import zipfile
import plot_fns
from emails import run_finished_notification
from app.get_effective_pathways import MutationTable, save_project_params_txt
from app.naming_rules import get_user_folder, get_project_folder, \
    get_params_path

_cleanup_thread = None


def remove_oldies():
    """Delete uploads and project directories over age limit."""
    max_age_days = current_app.config['PROJ_MAX_AGE_DAYS']
    cutoff = datetime.utcnow() - timedelta(days=max_age_days)
    oldies = UserFile.query.join(User).join(User.roles).\
        filter(Role.name == 'general').\
        filter(UserFile.upload_time < cutoff).all()
    old_ids = [userfile.file_id for userfile in oldies]
    for upload in oldies:
        if upload.run_complete:
            save_project_params_txt(upload)
            zip_project(upload)
        delete_project_folder(upload)
        db.session.delete(upload)
        db.session.commit()
    if oldies:
        current_app.logger.debug("Deleted project(s): {}".format(old_ids))


def delete_old_anonymous_users():
    max_age_days = current_app.config['ANONYMOUS_MAX_AGE_DAYS']
    cutoff = datetime.utcnow() - timedelta(days=max_age_days)
    # old_users = User.query.join(Role.users).\
    user_upload_tuples = db.session.query(User, UserFile).join(UserFile).\
        join(Role, User.roles).filter(Role.name == 'anonymous').\
        filter(User.member_since < cutoff).all()
    old_users = {i[0] for i in user_upload_tuples}
    old_uploads = {i[1] for i in user_upload_tuples}
    for upload in old_uploads:
        db.session.delete(upload)
    for user in old_users:
        user_folder = get_user_folder(user.id)
        if os.path.exists(user_folder):
            shutil.rmtree(user_folder)
        current_app.logger.debug("Deleting user {}.".format(user.email))
        db.session.delete(user)
    db.session.commit()


def delete_project_folder(upload):
    proj_folder = get_project_folder(upload)
    shutil.rmtree(proj_folder)


def start_cleanup_thread():
    global _cleanup_thread
    if _cleanup_thread is None:
        current_app.logger.info("Starting cleanup thread...")
        _cleanup_thread = Thread(target=tidy_projects_loop,
                                 args=[current_app._get_current_object()])
        _cleanup_thread.start()


def tidy_projects_loop(app, sleep_secs=None):
    if sleep_secs is None:
        sleep_secs = app.config['CLEANUP_INTERVAL']
    while True:
        time.sleep(sleep_secs)
        with app.app_context():
            remove_oldies()
            delete_old_anonymous_users()


def force_all_stopped_status():
    stop_list = UserFile.query.filter_by(run_complete=0).all()
    if stop_list:
        for upload_obj in stop_list:
            upload_id = upload_obj.file_id
            upload_obj.run_complete = None
            db.session.add(upload_obj)
            run_finished_notification(upload_id)
        db.session.commit()


def save_local_html(upload_obj, page_type='flat'):
    """Save local html file for {flat, scatter, tree} project results."""
    if page_type not in {'flat', 'scatter', 'tree'}:
        raise ValueError("page_type must be flat, scatter or tree.")

    if page_type == 'flat':
        output = get_flat_results(upload_obj)
    elif page_type == 'scatter':
        output = get_scatter_results(upload_obj)
    else:
        output = get_tree_results(upload_obj)

    proj_name = upload_obj.get_local_filename()
    proj_folder = get_project_folder(upload_obj)
    html_path = os.path.join(proj_folder, proj_name + '_' + page_type + '.html')
    with open(html_path, 'w') as out:
        out.write(output)


def get_flat_results(upload_obj):
    proj_names = {int(upload_obj.file_id): upload_obj.get_local_filename()}
    return render_template('pway/show_pathways_template.html',
                           proj=upload_obj,
                           projects=[upload_obj],
                           user_id=upload_obj.user_id,
                           show_proj=int(upload_obj.file_id),
                           proj_names=proj_names,
                           include_genes=None,
                           is_archive=True)


def get_scatter_results(upload_obj):
    """Render scatter results for specified project."""
    # get dictionary with js_name, js_inds, plot_inds, has_cnv, script, div
    scatter_dict = plot_fns.get_scatter_dict(upload_obj)
    return render_template('pway/scatter.html', current_proj=upload_obj,
                           projects=[upload_obj], user_id=upload_obj.user_id,
                           include_genes=None, is_archive=True, proj=upload_obj,
                           **scatter_dict)


def get_tree_results(upload_obj):
    """Render tree results for specified project."""
    proj_names = {int(upload_obj.file_id): upload_obj.get_local_filename()}
    names_odict = plot_fns.get_tree_data(upload_obj)
    return render_template('pway/tree.html', current_proj=upload_obj,
                           projects=[upload_obj], user_id=upload_obj.user_id,
                           proj_names=proj_names, names_odict=names_odict,
                           proj=upload_obj, is_archive=True)


def zip_project(upload_obj):
    """Zip project and save in user folder. Return path to zip file."""
    # generate params file if not already present
    params_path = get_params_path(upload_obj)
    if not os.path.exists(params_path):
        save_project_params_txt(upload_obj)
    # CREATE STATIC PAGES
    save_local_html(upload_obj, page_type='flat')
    save_local_html(upload_obj, page_type='scatter')
    save_local_html(upload_obj, page_type='tree')

    proj_name = upload_obj.get_local_filename()
    user_folder = get_user_folder(upload_obj.user_id)
    zip_path = os.path.join(user_folder, proj_name + '.zip')
    if os.path.exists(zip_path):
        os.remove(zip_path)
    zipf = zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED)
    proj_path = get_project_folder(upload_obj)
    for root, dirs, files in os.walk(proj_path):
        for this_file in files:
            # ignore raw .txt files in root of project directory
            if root is proj_path and \
                    (this_file.endswith('score_names.txt') or
                     this_file.endswith('scores.txt') or
                     this_file.endswith('.js') or
                     this_file.startswith('pathways')):
                continue
            zipf.write(os.path.join(root, this_file),
                       os.path.relpath(os.path.join(root, this_file),
                                       os.path.join(proj_path, '..')))
    zipf.close()
    return zip_path


def create_table(upload_obj):
    """Create table using uploaded data file."""
    proj_folder = get_project_folder(upload_obj)
    file_path = os.path.join(proj_folder, upload_obj.get_local_filename())
    table_name = upload_obj.get_table_name()
    table = MutationTable(table_name, file_path)
    if not table.loaded:
        raise Exception("Failed to load table {!r}.".format(table_name))
