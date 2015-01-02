from flask import render_template, flash, redirect, url_for, abort,\
    request, current_app
from . import pway, FileTester
from .forms import UploadForm
from ..models import UserFile
from .. import db
from flask_login import login_required, current_user
from werkzeug.utils import secure_filename
import os
from ..get_effective_pathways import run_analysis




@pway.route('/')
@login_required
def index():

    upload_list = UserFile.query.filter_by(user_id=current_user.id).all()
    return render_template('pway/index2.html', projects=upload_list)


@pway.route('/demo')
@login_required
def demo():
    return render_template('pway/show_pathways.html')

@pway.route('/results')
@login_required
def results():
    show_proj = request.args.get('proj', None)
    upload_list = UserFile.query.filter_by(user_id=current_user.id).filter_by(run_complete=True).all()
    return render_template('pway/show_pathways_template.html',
                           projects=upload_list, user_id=current_user.id,
                           show_proj=show_proj)


@pway.route('/upload', methods=('GET', 'POST'))
@login_required
def upload():
    """http://flask.pocoo.org/docs/0.10/patterns/fileuploads/"""
    upload_folder = current_app.config['UPLOAD_FOLDER']
    form = UploadForm()
    if form.validate_on_submit():
        # upload = Upload(uploader=current_user)
        mut_filename = secure_filename(form.mut_file.data.filename)

        # only accept this file if user has no running jobs.
        # will create folder for job: <run_id>
        user_upload = UserFile(filename=mut_filename, user_id=current_user.id)
        form.to_model(user_upload)
        db.session.add(user_upload)
        db.session.commit()
        user_folder = os.path.join(upload_folder, str(current_user.id))
        proj_folder = os.path.join(user_folder, str(user_upload.file_id))
        if not os.path.exists(user_folder):
            os.mkdir(user_folder)
        os.mkdir(proj_folder)
        file_path = os.path.join(proj_folder,
                                 user_upload.get_local_filename())
        form.mut_file.data.save(file_path)

        # FILE SHOULD BE IN UPLOADS FOLDER NOW.
        # RESPOND BASED ON FILE VALIDITY:
        #     VALID: START RUN, FLASH 'EXPECT EMAIL' on STATUS PAGE
        #     INVALID: FLASH 'BAD FILE', LIST ERRORS ON UPLOAD PAGE
        file_tester = FileTester(file_path)

        if file_tester.good_headers:
            if file_tester.good_headers and file_tester.data_present:
                user_upload.is_valid = True
                flash('File accepted and validated. Analysis in progress.')
                # want analysis to run asynchronously!
                # http://blog.miguelgrinberg.com/post/the-flask-mega-tutorial-part-xi-email-support
                run_analysis(proj_folder, file_path, user_upload)
                return redirect(url_for('.index'))
            else:  # good headers but no data
                flash('Your file seems to be missing data. Please try again.')
        else:  # bad headers
            flash("Your headers don't look right. Expected headers are:\n{}"
                .format('\t'.join(FileTester.want_headers)))
        db.session.add(user_upload)
        db.session.commit()
    return render_template('pway/upload.html', form=form)


"""

orig defaults:
 - output files saved to cwd
 - txt files: pway_pvalues_<PROJ_NAME>_<SUFFIX>(_pretty).txt
 - directories:
    - pathways_svg
    - matrix_txt
        - matrix_svg

UPLOADS FOLDER
==============

F 1_testfile.txt
F 2_melanoma_mk2014.txt
...
D output_1/
D output_2/
...


PROJECT RUN FOLDER
==================

 - create folder for project if text file is valid and user allowance not full

D pathways_svg/
D matrix_svg/
D matrix_txt/
F FILE_ID.txt
F FILE_ID_detail.txt
F FILE_ID.js

 - txt file name/path (default pathways_pvalues_PROJSTR) determined by
  GenericPathwayFileProcessor

root_svg_target = "pathways_svg/"
root_svg_matrix = "matrix_txt/matrix_svg/"
EXPECTED BY HTML: root_svg_matrix2 = "matrix_svg/"
proj_svg = 'proj_' + project_str + '/'

"""