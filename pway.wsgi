import os
import sys
import site
from glob import glob

abspath = os.path.dirname(__file__)  # /path/to/current/directory
os.chdir(abspath)
site.addsitedir(abspath)

from manage import app as application
application.debug = False

if __name__ == '__main__':
    application.debug = False
    application.run()
