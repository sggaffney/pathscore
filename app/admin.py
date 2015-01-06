from datetime import datetime, timedelta
from models import UserFile
from flask import current_app
from . import db
import os
import shutil
from threading import Thread
import time


_cleanup_thread = None

def get_user_folder(user_id):
    upload_folder = current_app.config['UPLOAD_FOLDER']
    return os.path.join(upload_folder, str(user_id))


def get_project_folder(upload_obj):
    user_folder = get_user_folder(upload_obj.user_id)
    return os.path.join(user_folder, str(upload_obj.file_id))


def remove_oldies():
    """Delete uploads and project directories over 3 weeks old."""
    cutoff = datetime.now() - timedelta(days=4, hours=0)
    oldies = UserFile.query.filter(UserFile.upload_time < cutoff).all()
    for upload in oldies:
        proj_folder = get_project_folder(upload)
        shutil.rmtree(proj_folder)
        db.session.delete(upload)
    db.session.commit()


def start_cleanup_thread():
    global _cleanup_thread
    if _cleanup_thread is None:
        print("Starting cleanup thread...")
        _cleanup_thread = Thread(target=tidy_projects_loop,
                                 args=[current_app._get_current_object()])
        _cleanup_thread.start()


def tidy_projects_loop(app):
    while True:
        time.sleep(app.config['PROJ_DELETE_INTERVAL'])
        with app.app_context():
            remove_oldies()