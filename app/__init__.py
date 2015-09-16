from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from flask_security import Security, SQLAlchemyUserDatastore  # for email verif.
from flask_login import LoginManager
from flask_bootstrap import Bootstrap
from flask_mail import Mail
from flask.ext.security.signals import user_registered
from config import config
import logging

bootstrap = Bootstrap()
db = SQLAlchemy()
mail = Mail()

login_manager = LoginManager()
# login_manager.login_view = 'auth.login'

from .models import User, Role

dbvars = dict()  # set in create_app. used in mysqldb queries.

def create_app(config_name):
    app = Flask(__name__)
    # development, testing, production, default
    app.config.from_object(config[config_name])

    logging_level = getattr(logging, app.config['LOGGING_LEVEL'],
                            'DEBUG')
    log_str = '[%(levelname)s] (%(threadName)-10s) %(message)s'
    logging.basicConfig(level=logging_level, format=log_str)

    if not app.debug:
        from logging.handlers import RotatingFileHandler
        # rotating 1MB log, with up to 10 backups
        file_handler = RotatingFileHandler(app.config['LOG_PATH'],
                                           'a', 1 * 1024 * 1024, 10)
        log_str = '%(asctime)s [%(levelname)s]: %(message)s [in %(pathname)s:%(lineno)d]'
        file_handler.setFormatter(logging.Formatter(log_str))
        app.logger.setLevel(logging.INFO)
        file_handler.setLevel(logging.INFO)
        app.logger.addHandler(file_handler)
        app.logger.info('Starting app.')

    global dbvars
    dbvars = dict(host=app.config['SGG_DB_HOST'],
                  db=app.config['SGG_DB_NAME'],
                  read_default_file=app.config['SGG_DB_CNF'])

    bootstrap.init_app(app)
    db.init_app(app)
    mail.init_app(app)

    from .pway import pway as pway_blueprint
    app.register_blueprint(pway_blueprint)

    from .auth import auth as auth_blueprint
    app.register_blueprint(auth_blueprint)

    from .momentjs import momentjs
    app.jinja_env.globals['momentjs'] = momentjs

    # Setup Flask-Security

    user_datastore = SQLAlchemyUserDatastore(db, User, Role)
    security = Security(app, user_datastore)
                        # login_form='auth.login')
                        # , confirm_register_form=None,
                        # register_form=None, forgot_password_form=None,
                        # reset_password_form=None,
                        # change_password_form=None, send_confirmation_form=None,
                        # passwordless_login_form=None)

    # security.init_app(app)

    # # Create a user to test with
    # @app.before_first_request
    # def create_user():
    #     db.create_all()
    #     user_datastore.create_user(email='sggaffney@gmail.com',
    #                                password='password')
    #     db.session.commit()

    from app.admin import start_cleanup_thread, force_all_stopped_status

    @app.before_first_request
    def before_first_request():
        start_cleanup_thread()
        force_all_stopped_status()

    @user_registered.connect_via(app)
    def user_registered_sighandler(app, user, confirm_token):
        default_role = user_datastore.find_role("general")
        user_datastore.add_role_to_user(user, default_role)
        db.session.commit()

    return app








