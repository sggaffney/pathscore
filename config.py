import os
import redis
from kombu import Exchange, Queue

basedir = os.path.abspath(os.path.dirname(__file__))

os.environ["PATH"] += os.pathsep + os.environ.get('SVGO_PARENT', '/usr/local/bin')

class Config:
    SECRET_KEY = os.environ.get('SECRET_KEY')

    MAIL_SERVER = os.environ.get('MAIL_SERVER')
    MAIL_PORT = int(os.environ.get('MAIL_PORT', 587))
    MAIL_USE_TLS = os.environ.get('MAIL_USE_TLS') == 'True'
    MAIL_USE_SSL = os.environ.get('MAIL_USE_SSL') == 'True'
    MAIL_USERNAME = os.environ.get('MAIL_USERNAME')
    MAIL_PASSWORD = os.environ.get('MAIL_PASSWORD')
    MAIL_DEFAULT_SENDER = os.environ.get('MAIL_DEFAULT_SENDER', '')
    MAIL_SENDER = os.environ.get('MAIL_SENDER')
    # MAIL_FLUSH_INTERVAL = 3600  # one hour
    MAIL_ERROR_RECIPIENT = os.environ.get('MAIL_ERROR_RECIPIENT')

    SECURITY_EMAIL_SENDER = os.environ.get('SECURITY_EMAIL_SENDER')
    SECURITY_PASSWORD_HASH = os.environ.get('SECURITY_PASSWORD_HASH')
    SECURITY_PASSWORD_SALT = os.environ.get('SECURITY_PASSWORD_SALT')
    SECURITY_CONFIRMABLE = os.environ.get('SECURITY_CONFIRMABLE') == 'True'
    SECURITY_REGISTERABLE = os.environ.get('SECURITY_REGISTERABLE') == 'True'
    SECURITY_CHANGEABLE = os.environ.get('SECURITY_CHANGEABLE') == 'True'
    SECURITY_EMAIL_SUBJECT_REGISTER = '[pathscore] Please confirm your account.'
    SECURITY_EMAIL_SUBJECT_CONFIRM = '[pathscore] Please confirm your account.'
    SECURITY_MSG_INVALID_PASSWORD = ("Bad username or password", "error")
    SECURITY_MSG_PASSWORD_NOT_PROVIDED = ("Bad username or password", "error")
    SECURITY_MSG_USER_DOES_NOT_EXIST = ("Bad username or password", "error")
    SECURITY_MSG_LOGIN = ('Please log in or select Register or Upload above. '
                          'You do not need an account to upload a project.',
                          'info')

    DATA_ROOT = os.environ.get('DATA_ROOT')
    TEMP_FOLDER = os.environ.get('TEMP_FOLDER')
    LOG_PATH = os.environ.get('LOG_PATH')

    MATRIX_BOX_PX = os.environ.get('MATRIX_BOX_PX', 25)
    MATRIX_HPAD_PX = os.environ.get('MATRIX_HPAD_PX', 220)
    MATRIX_TXT_SIZE = os.environ.get('MATRIX_TXT_SIZE', 10)

    DB_CNF = os.environ.get('MYSQLDB_CNF')
    DB_HOST = os.environ.get('MYSQLDB_HOST')
    # recycle connection before mysql default 8hr wait timeout
    SQLALCHEMY_POOL_RECYCLE = int(os.environ.get('SQLALCHEMY_POOL_RECYCLE', 3600))
    SQLALCHEMY_TRACK_MODIFICATIONS = False

    SERVER_NAME = os.environ.get('SERVER_NAME')

    CLEANUP_INTERVAL = 21600  # 6 hours
    PROJ_MAX_AGE_DAYS = 31
    ANONYMOUS_MAX_AGE_DAYS = 3

    LOGGING_LEVEL = 'DEBUG'

    USE_TOKEN_AUTH = True

    CELERY_BROKER_URL = 'redis://localhost:6379/0'
    CELERY_RESULT_BACKEND = 'redis://localhost:6379/0'
    CELERY_RESULT_DB_SHORT_LIVED_SESSIONS = True
    CELERY_DEFAULT_QUEUE = 'default'
    CELERY_CREATE_MISSING_QUEUES = True
    CELERY_ROUTES = {
        'app.emails.run_finished_notification_async':
            {'queue': 'mail'},
        'app.get_effective_pathways.run_analysis_async':
            {'queue': 'analysis'},
    }

    # enable rate limits only if redis is running
    try:
        r = redis.Redis()
        r.ping()
        USE_RATE_LIMITS = True
    except redis.ConnectionError:
        USE_RATE_LIMITS = False


class DevelopmentConfig(Config):
    DEBUG = True
    SECRET_KEY = os.environ.get('SECRET_KEY') or 't0p s3cr3t'
    SQLALCHEMY_DATABASE_URI = os.environ.get('DEV_DATABASE_URL') or \
        'sqlite:///' + os.path.join(basedir, 'data-dev.sqlite')
    MAIL_FLUSH_INTERVAL = 60  # one minute
    DB_NAME = os.environ.get('MYSQLDB_DB_DEV')
    SQLALCHEMY_ECHO = os.environ.get('SQLALCHEMY_ECHO') == 'True'

class TestingConfig(Config):
    TESTING = True
    SECRET_KEY = 'secret'
    SQLALCHEMY_DATABASE_URI = os.environ.get('TEST_DATABASE_URL') or \
        'sqlite:///' + os.path.join(basedir, 'data-test.sqlite')
    DB_NAME = os.environ.get('MYSQLDB_DB_TEST')


class ProductionConfig(Config):
    SQLALCHEMY_DATABASE_URI = os.environ.get('DATABASE_URL') or \
        'sqlite:///' + os.path.join(basedir, 'data.sqlite')
    DB_NAME = os.environ.get('MYSQLDB_DB')


config = {
    'development': DevelopmentConfig,
    'testing': TestingConfig,
    'production': ProductionConfig,

    'default': DevelopmentConfig
}


