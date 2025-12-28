import os
import logging
from dotenv import load_dotenv, find_dotenv
from typing import Union


def _get_env_path():
    """Use .env in app project directory, fallback to find_dotenv.

    Returns None if no .env file is found (e.g., when running in Docker
    with environment variables passed via docker-compose).
    """
    basedir = os.path.abspath(os.path.dirname(__file__))
    default_path = os.path.join(basedir, '.env')
    if os.path.exists(default_path):
        return default_path
    alt_path = find_dotenv()
    if alt_path:
        return alt_path
    return None  # No .env file found - use environment variables directly


def _create_logger(log_level: Union[str, None] = None):
    if log_level is None:
        log_level = os.environ.get('LOG_LEVEL', 'INFO').upper()
    log_level_int = getattr(logging, log_level)

    ch = logging.StreamHandler()
    ch.setLevel(log_level_int)
    formatter = logging.Formatter(
        fmt='%(asctime)s - %(message)s',
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    ch.setFormatter(formatter)

    _logger = logging.getLogger(__name__)
    _logger.addHandler(ch)
    _logger.setLevel(log_level_int)

    logging.getLogger('pway_app.app').setLevel(log_level_int)
    return _logger


env_path = _get_env_path()
if env_path:
    print(f"Loading .env from {env_path}")
    load_dotenv(env_path, override=False)
else:
    print("No .env file found - using environment variables directly")

logger = _create_logger()


from app import create_app

app = create_app(os.getenv('FLASK_CONFIG') or 'default')
