from datetime import datetime, timedelta
from models import UserFile, User, Role
from flask import current_app
from . import db
import os
import shutil
from threading import Thread
import time
import zipfile
from emails import run_finished_notification
from app.get_effective_pathways import MutationTable
from app.naming_rules import get_user_folder, get_project_folder


_cleanup_thread = None


def remove_oldies():
    """Delete uploads and project directories over age limit."""
    max_age_days = current_app.config['PROJ_MAX_AGE_DAYS']
    cutoff = datetime.now() - timedelta(days=max_age_days)
    oldies = UserFile.query.join(User).join(User.roles).\
        filter(Role.name == 'general').\
        filter(UserFile.upload_time < cutoff).all()
    for upload in oldies:
        if upload.run_complete:
            save_project_params_txt(upload)
            zip_project(upload)
        delete_project_folder(upload)
        db.session.delete(upload)
        db.session.commit()


def delete_old_anonymous_users():
    max_age_days = current_app.config['ANONYMOUS_MAX_AGE_DAYS']
    cutoff = datetime.now() - timedelta(days=max_age_days)
    # old_users = User.query.join(Role.users).\
    user_upload_tuples = db.session.query(User, UserFile).join(UserFile).\
        join(Role, User.roles).filter(Role.name=='vip').\
        filter(User.member_since < cutoff).all()
    old_users = {i[0] for i in user_upload_tuples}
    old_uploads = {i[1] for i in user_upload_tuples}
    for upload in old_uploads:
        db.session.delete(upload)
    for user in old_users:
        user_folder = get_user_folder(user.id)
        if os.path.exists(user_folder):
            shutil.rmtree(user_folder)
        db.session.delete(user)
    db.session.commit()



def delete_project_folder(upload):
    proj_folder = get_project_folder(upload)
    shutil.rmtree(proj_folder)


def start_cleanup_thread():
    global _cleanup_thread
    if _cleanup_thread is None:
        print("Starting cleanup thread...")
        _cleanup_thread = Thread(target=tidy_projects_loop,
                                 args=[current_app._get_current_object()])
        _cleanup_thread.start()


def tidy_projects_loop(app):
    while True:
        time.sleep(app.config['CLEANUP_INTERVAL'])
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



# def zipdir(path, zip):
#     """Called by zip_project."""
#     for root, dirs, files in os.walk(path):
#         for file in files:
#             zip.write(os.path.join(root, file))


def zip_project(upload_obj):
    # TODO: only zip if file not there already
    proj_name = upload_obj.get_local_filename()
    user_folder = get_user_folder(upload_obj.user_id)
    zip_path = os.path.join(user_folder, proj_name + '.zip')
    zipf = zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED)
    proj_path = get_project_folder(upload_obj)
    for root, dirs, files in os.walk(proj_path):
        for file in files:
            # ignore raw .txt files in root of project directory
            if root is proj_path and file.endswith('.txt'):
                continue
            zipf.write(os.path.join(root, file),
                       os.path.relpath(os.path.join(root, file),
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
