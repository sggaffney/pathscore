from threading import Thread
from functools import wraps
from flask import request, redirect

def no_ssl(fn):
    @wraps(fn)
    def decorated_view(*args, **kwargs):
        if request.is_secure:
            return redirect(request.url.replace("https://", "http://"))
        else:
            return fn(*args, **kwargs)
        return fn(*args, **kwargs)
    return decorated_view

def async(fn):
    @wraps(fn)
    def wrapper(*args, **kwargs):
        thr = Thread(target=fn, args=args, kwargs=kwargs)
        thr.start()
    return wrapper
