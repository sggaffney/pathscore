from threading import Thread
from functools import wraps
from datetime import datetime, timedelta
from flask import request, redirect, g
from flask_login import current_user
from models import UserFile
from misc import get_wait_time_string
from errors import LimitError


def no_ssl(fn):
    @wraps(fn)
    def decorated_view(*args, **kwargs):
        if request.is_secure:
            return redirect(request.url.replace("https://", "http://"))
        else:
            return fn(*args, **kwargs)
    return decorated_view


def ssl_required(fn):
    @wraps(fn)
    def decorated_view(*args, **kwargs):
        if request.is_secure:
            return fn(*args, **kwargs)
        else:
            return redirect(request.url.replace("http://", "https://"))
        return fn(*args, **kwargs)

    return decorated_view


def async(fn):
    @wraps(fn)
    def wrapper(*args, **kwargs):
        thr = Thread(target=fn, args=args, kwargs=kwargs)
        thr.start()
    return wrapper


def limit_user_uploads(f):
    """This decorator implements upload limiting."""
    @wraps(f)
    def wrapped(*args, **kwargs):
        # CHECK IF RUNNING PROJECT COUNT IS WITHIN USER LIMITS
        if current_user.is_authenticated:
            incomplete = UserFile.query.filter_by(user_id=current_user.id)\
                .filter_by(run_complete=0).all()
            role_names = [r.name for r in current_user.roles]
            if 'vip' not in role_names and incomplete:
                raise LimitError("Sorry, you must wait until your currently "
                                 "running projects have finished.")

            # ENFORCE WEEKLY LIMIT AND DISPLAY INFO MESSAGE
            week_ago = datetime.utcnow() - timedelta(days=7)
            week_complete = UserFile.query.filter_by(user_id=current_user.id)\
                .filter_by(run_complete=True)\
                .filter(UserFile.upload_time > week_ago).all()
            n_week = len(week_complete)
            n_week_max = int(max([r.uploads_pw for r in current_user.roles]))
            # if len(week_complete>9):
            if n_week >= n_week_max:
                first_time = min([p.upload_time for p in week_complete])
                wait_time = first_time + timedelta(days=7) - datetime.utcnow()
                wait_str = get_wait_time_string(wait_time)
                raise LimitError("Sorry, you've used your allotted runs for "
                                 "this week. You can try again in {}."
                                 .format(wait_str))
            # ADD TO REQUEST GLOBAL, FOR USER FEEDBACK
            g.n_week = n_week
            g.n_week_max = n_week_max
            g.incomplete = incomplete

        # let the request through
        return f(*args, **kwargs)
    return wrapped
