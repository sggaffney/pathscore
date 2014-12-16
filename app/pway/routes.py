from flask import render_template, flash, redirect, url_for, abort,\
    request, current_app
from . import pway, FileTester
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
        file_path = os.path.join(current_app.config['UPLOAD_FOLDER'],
                                 user_upload.get_local_filename())
        form.mut_file.data.save(file_path)

        # FILE SHOULD BE IN UPLOADS FOLDER NOW.
        # RESPOND BASED ON FILE VALIDITY:
        #     VALID: START RUN, FLASH 'EXPECT EMAIL' on STATUS PAGE
        #     INVALID: FLASH 'BAD FILE', LIST ERRORS ON UPLOAD PAGE
        file_tester = FileTester(file_path)

        if file_tester.good_headers:
            if file_tester.good_headers and file_tester.data_present:
                flash('File accepted and validated. Analysis in progress.')
                return redirect(url_for('.index'))
            else:  # good headers but no data
                flash('Your file seems to be missing data. Please try again.')
        else:  # bad headers
            flash("Your headers don't look right. Expected headers are:\n{}"
                .format('\t'.join(FileTester.want_headers)))
    return render_template('pway/upload.html', form=form)
