from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from flask_security import Security, SQLAlchemyUserDatastore  # for email verif.
from flask_login import LoginManager
from flask_bootstrap import Bootstrap
from flask_mail import Mail
from flask.ext.security.signals import user_registered
from celery import Celery
from config import config, Config
import logging

bootstrap = Bootstrap()
db = SQLAlchemy()
mail = Mail()
celery = Celery(__name__, broker=Config.CELERY_BROKER_URL)

login_manager = LoginManager()
# login_manager.login_view = 'auth.login'

from .models import User, Role

dbvars = dict()  # set in create_app. used in mysqldb queries.


def create_app(config_name):
    app = Flask(__name__)
    # development, testing, production, default
    app.config.from_object(config[config_name])

    logging_level = getattr(logging, app.config['LOGGING_LEVEL'], logging.DEBUG)
    log_str = ('%(asctime)s [%(levelname)s]: %(message)s '
               '[(%(threadName)-10s) in %(pathname)s:%(lineno)d]')
    if not app.debug:
        # FILE LOG [info]. rotating 1MB log, with up to 10 backups
        from logging.handlers import RotatingFileHandler
        file_handler = RotatingFileHandler(app.config['LOG_PATH'],
                                           'a', 1 * 1024 * 1024, 10)
        file_handler.setFormatter(logging.Formatter(log_str))
        file_handler.setLevel(logging_level)
        app.logger.addHandler(file_handler)
        # EMAIL LOG [error].
        from logging.handlers import SMTPHandler
        mail_handler = SMTPHandler(app.config['MAIL_SERVER'],
                                   app.config['MAIL_SENDER'],
                                   app.config['MAIL_ERROR_RECIPIENT'],
                                   '[pathscore] Error report')
        mail_handler.setFormatter(logging.Formatter(
            "Msg type:  %(levelname)s\n"
            "Location:  %(pathname)s:%(lineno)d\n"
            "Module:    %(module)s\n"
            "Function:  %(funcName)s\n"
            "Time:      %(asctime)s\n"
            "\n"
            "Message:\n"
            "\n"
            "%(message)s\n"))
        mail_handler.setLevel(logging.ERROR)
        app.logger.addHandler(mail_handler)
    else:
        logging.basicConfig(format=log_str, level=logging_level)
    app.logger.setLevel(logging.INFO)
    # LOG STARTUP
    app.logger.info('Starting app')

    global dbvars
    dbvars = dict(host=app.config['DB_HOST'],
                  db=app.config['DB_NAME'],
                  read_default_file=app.config['DB_CNF'])

    bootstrap.init_app(app)
    db.init_app(app)
    mail.init_app(app)
    celery.conf.update(app.config)

    from .pway import pway as pway_blueprint
    app.register_blueprint(pway_blueprint)

    from .auth import auth as auth_blueprint
    app.register_blueprint(auth_blueprint)

    from .api import api as api_blueprint
    app.register_blueprint(api_blueprint, url_prefix='/api')

    if app.config['USE_TOKEN_AUTH']:
        from api.token import token as token_blueprint
        app.register_blueprint(token_blueprint, url_prefix='/auth')

    from .momentjs import momentjs
    app.jinja_env.globals['momentjs'] = momentjs

    # Setup Flask-Security

    user_datastore = SQLAlchemyUserDatastore(db, User, Role)
    Security(app, user_datastore)
             # login_form='auth.login', confirm_register_form=None,
             # register_form=None, forgot_password_form=None,
             # reset_password_form=None, change_password_form=None,
             # send_confirmation_form=None, passwordless_login_form=None)

    # # Create a user to test with
    # @app.before_first_request
    # def create_user():
    #     db.create_all()
    #     user_datastore.create_user(email='fake@yale.edu',
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
        app.logger.info('Registered user {}'.format(user.email))
        db.session.commit()

    @app.context_processor
    def inject_css():
        return dict(get_common_css=get_common_css)

    return app


def get_common_css(is_archive=False):
    """Get text of common_js or script tag with link."""
    if not is_archive:
        css_text = '<link rel="stylesheet" type="text/css" href="' \
                   + url_for('static', filename='styles.css') + '">'
    else:
        css_path = os.path.join(current_app.static_folder,
                                'styles.css')
        with open(css_path, 'r') as infile:
            css_text = infile.read()
        css_text = "<style>\n" + css_text + "\n</style>"
    return css_text
