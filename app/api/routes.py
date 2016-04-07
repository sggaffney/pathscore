import os
from flask import current_app, g, send_file  # redirect, url_for, abort,
from flask_login import current_user, login_user
from flask import request
from . import api
from app import db
from ..decorators import limit_user_uploads  # ssl_required
from decorators import json, etag, collection
from ..errors import ValidationError
from ..uploads import MutationFile
from ..models import UserFile, create_anonymous_user, initialize_project
from ..admin import delete_project_folder
from ..get_effective_pathways import run_analysis
from ..admin import zip_project
from auth import auth, auth_optional


@api.route('/', methods=['POST'])
# @ssl_required
@auth.login_required
@json
def test():
    a = request.get_json()
    a.update({'hola': 'amigo'})
    return a, 201, {'Location': 'some_link.html'}


@api.route('/archives/<int:proj>', methods=['GET'])
@auth.login_required
def archive(proj):
    upload_obj = UserFile.query.\
        filter_by(user_id=g.user.id, file_id=proj).\
        first_or_404()
    zip_path = zip_project(upload_obj)
    filename = os.path.basename(zip_path)
    return send_file(zip_path, mimetype='application/zip',
                     as_attachment=True, attachment_filename=filename)


@api.route('/projects/', methods=['GET'])
@etag
@auth.login_required
@json
@collection(UserFile, name='projects')
def get_user_projects():
    return UserFile.query.filter_by(user_id=g.user.id)


@api.route('/projects/<int:file_id>', methods=['GET'])
@etag
@auth.login_required
@json
def get_project(file_id):
    return UserFile.query.filter_by(user_id=g.user.id, file_id=file_id).\
        first_or_404()


@api.route('/projects/<int:file_id>', methods=['DELETE'])
@etag
@auth.login_required
@json
def delete_project(file_id):
    project = UserFile.query.filter_by(user_id=g.user.id, file_id=file_id).\
        first_or_404()
    delete_project_folder(project)
    db.session.delete(project)
    db.session.commit()
    return {}


@api.route('/projects/', methods=['POST'])
@limit_user_uploads
# @ssl_required
# @auth_optional.login_required  <- IN BEFORE_REQUEST
@json
def upload():
    """http://flask.pocoo.org/docs/0.10/patterns/fileuploads/"""

    # VALIDATE FORM/FILE DATA
    user_upload = UserFile().import_data(request.form)
    filestore = request.files['mut_file']
    mut_filename = filestore.filename
    if not mut_filename.endswith('.txt') and not mut_filename.endswith('.tsv'):
        raise ValidationError("Use txt or tsv extension for mut_file.")
    mut_file = MutationFile(filestore)

    rv = {}
    # CREATE NEW USER IF UNAUTHENTICATED (HERE, FILE IS VALID)
    if not current_user.is_authenticated():
        # create guest user
        temp_user, temp_pswd = create_anonymous_user()
        rv['user_name'] = temp_user.email
        rv['user_password'] = temp_pswd
        rv['message'] = 'This temporary account will be deleted in {} days.'.\
            format(current_app.config['ANONYMOUS_MAX_AGE_DAYS'])
        login_user(temp_user, force=True, remember=True)

    # CREATE USERFILE OBJECT
    user_upload.filename = mut_filename
    user_upload.user_id = current_user.id
    out = initialize_project(user_upload=user_upload, mut_file=mut_file)
    user_upload, proj_folder, file_path = out

    success_msg = 'File accepted and validated. Analysis in progress.'
    rv['status'] = 'Success.'
    rv['message'] = ' '.join([success_msg, rv['message']]) if 'message' in rv \
        else success_msg

    # RUN ANALYSIS:
    run_analysis(proj_folder, file_path, user_upload.file_id)

    return rv, 201, {'Location': user_upload.get_url()}
