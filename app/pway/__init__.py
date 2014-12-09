from flask import Blueprint

pway = Blueprint('pway', __name__)

from . import routes

