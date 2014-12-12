from flask import render_template, flash, redirect, url_for, abort,\
    request, current_app
from . import pway
from .forms import UploadForm
from flask_login import login_required  # current_user
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
        form.mut_file.data.save(os.path.join(current_app.config['UPLOAD_FOLDER'],
                                             mut_filename))

        # form.to_model(upload)
        # db.session.add(upload)
        # db.session.commit()
        flash('The input files were added successfully.')
        return redirect(url_for('.index'))
    else:
        mut_filename = None
    return render_template('pway/upload.html', form=form, filename=mut_filename)

