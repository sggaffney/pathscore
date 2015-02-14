from flask import render_template, flash, redirect, url_for, abort,\
    request, current_app, send_file, send_from_directory
from . import pway, FileTester
from .forms import UploadForm
from ..models import UserFile
from .. import db
from flask_login import login_required, current_user
from werkzeug.utils import secure_filename
import os
from ..get_effective_pathways import run_analysis
from ..admin import get_project_folder, get_user_folder, zip_project, \
    delete_project_folder
from datetime import datetime, timedelta

@pway.route('/')
@login_required
def index():

    upload_list = UserFile.query.filter_by(user_id=current_user.id).all()
    return render_template('pway/index2.html', projects=upload_list)


@pway.route('/demo')
@login_required
def demo():
    return render_template('pway/show_pathways_demo.html')


@pway.route('/archive')
@login_required
def archive():
    proj_id = int(request.args.get('proj', 0))
    upload_obj = None
    if proj_id:
        upload_obj = UserFile.query.get(proj_id)
    if not upload_obj or current_user.id != upload_obj.user_id:
        abort(404)
    zip_path = zip_project(upload_obj)
    filename = os.path.basename(zip_path)
    return send_file(zip_path, mimetype='application/zip',
                     as_attachment=True, attachment_filename=filename)

@pway.route('/demofile')
@login_required
def demo_file():
    # return render_template('pway/show_pathways_template.html')
    return send_from_directory(current_app.config['UPLOAD_FOLDER'],
                               'skcm_ns_500.txt', as_attachment=True)


@pway.route('/results')
@login_required
def results():
    show_proj = request.args.get('proj', None)
    upload_list = UserFile.query.filter_by(user_id=current_user.id).filter_by(run_complete=True).all()
    proj_names = {int(i.file_id): i.get_local_filename() for i in upload_list}
    if not(upload_list):
        flash("No project results to show yet.", "info")
    return render_template('pway/show_pathways_template.html',
                           projects=upload_list, user_id=current_user.id,
                           show_proj=show_proj, proj_names=proj_names)


@pway.route('/upload', methods=('GET', 'POST'))
@login_required
def upload():
    """http://flask.pocoo.org/docs/0.10/patterns/fileuploads/"""
    # only accept this file if user has no running jobs.
    user_id = current_user.id
    incomplete = UserFile.query.filter_by(user_id=user_id)\
        .filter_by(run_complete=0).all()
    role_names = [r.name for r in current_user.roles]
    if 'townsend' not in role_names and incomplete:
        flash("Sorry, you must wait until your currently running projects "
              "have finished.", "danger")
        return index()

    # ENFORCE WEEKLY LIMIT AND DISPLAY INFO MESSAGE
    week_ago = datetime.now() - timedelta(days=7)
    week_complete = UserFile.query.filter_by(user_id=user_id)\
        .filter_by(run_complete=True).filter(UserFile.upload_time > week_ago)\
        .all()
    n_week = len(week_complete)
    n_week_max = int(max([r.uploads_pw for r in current_user.roles]))
    # if len(week_complete>9):
    if n_week >= n_week_max:
        first_time = min([p.upload_time for p in week_complete])
        wait_time = first_time + timedelta(days=7) - datetime.now()
        diff_str = '{:.1f}'.format(wait_time.seconds/(60*60.))

        flash("Sorry, you've used your allotted runs for this week. You can "
              "try again in {} hours.".format(diff_str), "danger")
        return index()

    form = UploadForm()
    if form.validate_on_submit():
        # upload = Upload(uploader=current_user)
        mut_filename = secure_filename(form.mut_file.data.filename)

        # will create folder for job: <run_id>
        user_upload = UserFile(filename=mut_filename, user_id=user_id)
        form.to_model(user_upload)
        db.session.add(user_upload)
        db.session.commit()

        user_folder = get_user_folder(user_id)
        proj_folder = get_project_folder(user_upload)
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
                user_upload.run_complete = False
                upload_id = user_upload.file_id
                db.session.add(user_upload)
                db.session.commit()
                flash('File accepted and validated. Analysis in progress.',
                      'success')
                # want analysis to run asynchronously!
                # http://blog.miguelgrinberg.com/post/the-flask-mega-tutorial-part-xi-email-support
                run_analysis(proj_folder, file_path, upload_id)
                return redirect(url_for('.index'))
            else:  # good headers but no data
                flash('Your file seems to be missing data. Please try again.',
                      'danger')
        else:  # bad headers
            flash("Your headers don't look right. Expected headers are:\n{}"
                  .format('\t'.join(FileTester.want_headers)), 'danger')
            delete_project_folder(user_upload)
            db.session.delete(user_upload)
            db.session.commit()
    else:
        message = "You've run {} projects in the last week. ".format(n_week)
        message += "Your weekly limit is {}.".format(n_week_max)
        if incomplete:
            message += "\nYou have {} jobs still running."\
                .format(len(incomplete))
        flash(message, "info")
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