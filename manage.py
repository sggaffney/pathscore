#!/usr/bin/env python
import os
if os.path.exists('.env'):
    print('Importing environment from .env...')
    for line in open('.env'):
        var = line.strip().split('=')
        if len(var) == 2:
            os.environ[var[0]] = var[1]

from app import create_app
from flask_script import Manager
from app import db
from app.models import User, Role
from flask.ext.migrate import Migrate, MigrateCommand

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


if __name__ == '__main__':
    manager.run()

