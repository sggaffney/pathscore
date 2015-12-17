from flask import Blueprint, g, flash, url_for
from app.get_effective_pathways import path_info_dict
from ..errors import LimitError


pway = Blueprint('pway', __name__)


@pway.context_processor
def inject_n_pathways():
    return dict(n_pathways=len(path_info_dict))


@pway.before_request
def get_n_pathways():
    g.n_pathways = len(path_info_dict)


@pway.errorhandler(LimitError)
def validation_error(e):
    flash(str(e), "danger")
    return url_for('.index')


from . import routes
