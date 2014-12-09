from flask import render_template, flash, redirect, url_for, abort,\
    request, current_app
from . import pway
from flask_login import login_required, current_user

@pway.route('/')
@login_required
def index():
    # page = request.args.get('page', 1, type=int)
    # pagination = Talk.query.order_by(Talk.date.desc()).paginate(
    #     page, per_page=current_app.config['TALKS_PER_PAGE'],
    #     error_out=False)
    # talk_list = pagination.items
    return render_template('pway/show_pathways.html')
