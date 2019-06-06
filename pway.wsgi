import os
import sys
import site
from glob import glob

from dotenv import load_dotenv


load_dotenv()
abspath = os.path.dirname(__file__)  # /path/to/current/directory
os.chdir(abspath)
site.addsitedir(abspath)
sys.path.insert(1, os.getenv("SITE_PACKAGES_DIR"))

from manage import app as application
application.debug = False

if __name__ == '__main__':
    application.debug = False
    application.run()
