"""Assumes venv in ~/pway_venv and app in Downloads/pway_app."""

import os
import sys
import site

abspath = os.path.dirname(__file__)  # /Users/sgg/Downloads/pway_app
user_dir = abspath[0:abspath.index('/', 7)+1]
os.chdir(abspath)
site.addsitedir(abspath)
site.addsitedir(os.path.join(user_dir,'pway_venv/lib/python2.7/site-packages'))
#sys.path.insert(0, os.path.join(user_dir, 'pway_venv/lib/python2.7/site-packages'))
#sys.path.insert(0, '/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7')  # for building cython
#sys.path.insert(0,'/System/Library/Frameworks/Python.framework/Versions/2.7/Extras/lib/python/')
os.environ["PATH"] += os.pathsep + "/Applications/MATLAB_R2013a.app/bin"

from manage import app as application
application.debug = False

if __name__ == '__main__':
    application.debug = False
    application.run()

