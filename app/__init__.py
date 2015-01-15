from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from flask_security import Security, SQLAlchemyUserDatastore  # for email verif.
from flask_login import LoginManager
from flask_bootstrap import Bootstrap
from flask_mail import Mail
from flask.ext.security.signals import user_registered
from config import config


bootstrap = Bootstrap()
db = SQLAlchemy()
mail = Mail()

login_manager = LoginManager()
# login_manager.login_view = 'auth.login'

from .models import User, Role

def create_app(config_name):
    app = Flask(__name__)
    # development, testing, production, default
    app.config.from_object(config[config_name])

    app.config['SECURITY_CONFIRMABLE'] = True
    app.config['SECURITY_REGISTERABLE'] = True
    app.config['SECURITY_CHANGEABLE'] = True

    bootstrap.init_app(app)
    db.init_app(app)
    mail.init_app(app)

    from .pway import pway as pway_blueprint
    app.register_blueprint(pway_blueprint)

    from .auth import auth as auth_blueprint
    app.register_blueprint(auth_blueprint)

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








