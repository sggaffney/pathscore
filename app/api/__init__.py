from flask import Blueprint, g, url_for

from .errors import bad_request, not_found, not_allowed, too_many_requests
from ..errors import ValidationError, LimitError
from .auth import auth, auth_optional
from .decorators import json, rate_limit


api = Blueprint('api', __name__)


@api.errorhandler(ValidationError)
def validation_error(e):
    return bad_request(str(e))


@api.errorhandler(LimitError)
def limit_error(e):
    return too_many_requests(str(e))


@api.errorhandler(400)
def bad_request_error(e):
    return bad_request('invalid request')


@api.errorhandler(404)
@auth.login_required
def not_found_error(e):
    return not_found('item not found')


@api.errorhandler(405)
def method_not_allowed_error(e):
    return not_allowed()


@api.before_request
@auth_optional.login_required
@rate_limit(limit=5, period=15)
def before_request():
    pass


@api.after_request
def after_request(response):
    if hasattr(g, 'headers'):
        response.headers.extend(g.headers)
    return response

# do this last to avoid circular dependencies
from . import routes
