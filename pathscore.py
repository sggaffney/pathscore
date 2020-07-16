import os
import logging
from dotenv import load_dotenv, find_dotenv


logger = logging.getLogger(__name__)

# LOAD .env file from location specified by ENV_NAME environment variable
env_path = os.getenv('ENV_NAME', find_dotenv())
logger.info("Loading .env from %s", env_path)
load_dotenv(env_path, override=True)

from .app import create_app

app = create_app(os.getenv('FLASK_CONFIG') or 'default')
