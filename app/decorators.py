from threading import Thread
from functools import wraps
from flask import request, redirect

def no_ssl(fn):
    @wraps(fn)
    def decorated_view(*args, **kwargs):
        if request.is_secure:
            return redirect(request.url.replace("http://", "https://"))
        else:
            return fn(*args, **kwargs)
        return fn(*args, **kwargs)

    return decorated_view



def async(f):
    def wrapper(*args, **kwargs):
        thr = Thread(target=f, args=args, kwargs=kwargs)
        thr.start()
    return wrapper
