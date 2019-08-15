#!/usr/bin/env python
import os
import logging
from dotenv import load_dotenv, find_dotenv


logger = logging.getLogger(__name__)

# LOAD .env file from location specified by ENV_NAME environment variable
env_path = os.getenv('ENV_NAME', find_dotenv())
logger.info("Loading .env from %s", env_path)
load_dotenv(env_path, override=True)


import shutil
from flask_migrate import Migrate, MigrateCommand
from flask_script import Manager

from app import create_app, db, mds
from app.models import User, Role, UserFile
from app.naming_rules import get_project_folder


app = create_app(os.getenv('FLASK_CONFIG') or 'default')
migrate = Migrate(app, db)

manager = Manager(app)
manager.add_command('db', MigrateCommand)


@manager.command
def test():
    from subprocess import call
    call(['nosetests', '-v',
          '--with-coverage', '--cover-package=app', '--cover-branches',
          '--cover-erase', '--cover-html', '--cover-html-dir=cover'])


@manager.command
def adduser(email, role='general'):
    """Register a new user."""
    from getpass import getpass
    password = getpass()
    password2 = getpass(prompt='Confirm: ')
    if password != password2:
        import sys
        sys.exit('Error: passwords do not match.')
    role_obj = Role.query.filter_by(name=role).one()
    db.create_all()
    user = User(email=email, password_raw=password)
    user.roles.append(role_obj)
    # app.user_datastore.add_role_to_user(user, default_role)
    db.session.add(user)
    db.session.commit()
    print('User {0} was registered successfully.'.format(email))


@manager.option('-p', '--proj', dest='proj_str', default=None)
def delete_proj(proj_str):

    in_list = [i for i in proj_str.split(',')]
    proj_list = [int(i) for i in in_list if i.isdigit()]
    if len(in_list) != len(proj_list):
        import sys
        sys.exit("Provide project ids as comma-separated integers.")
    for proj in proj_list:
        upload = UserFile.query.get(proj)
        if upload:
            try:
                proj_folder = get_project_folder(upload)
                shutil.rmtree(proj_folder)
            except OSError as e:
                print(e)
            db.session.delete(upload)
            print("Project {} deleted successfully.".format(proj))
        else:
            print("Project {} not found.".format(proj))
    db.session.commit()


@manager.option('-p', '--proj', dest='proj_id', default=None)
def build_mds(proj_id):
    """Create MDS dataframe and points files."""
    mds.build_proj_mds(proj_id)


if __name__ == '__main__':
    manager.run()
