from flask import current_app, g
from flask_httpauth import HTTPBasicAuth
from flask_login import login_user
from ..models import User
from .errors import unauthorized


auth = HTTPBasicAuth()
auth_optional = HTTPBasicAuth()


@auth_optional.error_handler
@auth.error_handler
def unauthorized_error():
    return unauthorized()


@auth.verify_password
def verify_password(email_or_token, password):
    if current_app.config['USE_TOKEN_AUTH']:
        # token authentication
        g.user = User.validate_auth_token(email_or_token)
        if g.user is not None:
            login_user(g.user, force=True, remember=True)
            return True
    else:
        # username/password authentication
        g.user = User.query.filter_by(email=email_or_token).first()
        if g.user is not None and g.user.verify_password(password):
            login_user(g.user, force=True, remember=True)
            return True
    return False


@auth_optional.verify_password
def verify_optional_password(email_or_token='', password=''):
    if email_or_token == '':
        return True

    if current_app.config['USE_TOKEN_AUTH']:
        # token authentication
        g.user = User.validate_auth_token(email_or_token)
        if g.user is not None:
            login_user(g.user, force=True, remember=True)
            return True
    else:
        # username/password authentication
        g.user = User.query.filter_by(email=email_or_token).first()
        if g.user is not None and g.user.verify_password(password):
            login_user(g.user, force=True, remember=True)
            return True
    return False
