from flask import render_template, flash, redirect, url_for, abort,\
    request, current_app
from . import pway
from .forms import UploadForm
from ..models import UserFile
from .. import db
from flask_login import login_required, current_user
from werkzeug.utils import secure_filename
import os


@pway.route('/')
@login_required
def index():
    return render_template('pway/index2.html')


@pway.route('/demo')
@login_required
def demo():
    return render_template('pway/show_pathways.html')


@pway.route('/upload', methods=('GET', 'POST'))
@login_required
def upload():
    """http://flask.pocoo.org/docs/0.10/patterns/fileuploads/"""
    form = UploadForm()
    if form.validate_on_submit():
        # upload = Upload(uploader=current_user)
        mut_filename = secure_filename(form.mut_file.data.filename)

        # save mut_filename, time, size, user_id in uploads table
        # with 'is_valid' field, 'run_accepted', 'run_complete'
        # run_id primary key
        # only accept this file if user has no running jobs.
        # will create folder for job: <run_id>
        user_upload = UserFile(filename=mut_filename, user_id=current_user.id)
        form.to_model(user_upload)
        db.session.add(user_upload)
        db.session.commit()

        form.mut_file.data.save(
            os.path.join(current_app.config['UPLOAD_FOLDER'], mut_filename))

        # FILE SHOULD BE IN UPLOADS FOLDER NOW.
        # RESPOND BASED ON FILE VALIDITY:
        #     VALID: START RUN, FLASH 'EXPECT EMAIL' on STATUS PAGE
        #     INVALID: FLASH 'BAD FILE', LIST ERRORS ON UPLOAD PAGE


        flash('The input files were added successfully.')
        return redirect(url_for('.index'))
    return render_template('pway/upload.html', form=form)

