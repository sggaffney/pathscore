import os
basedir = os.path.abspath(os.path.dirname(__file__))


# 	private final String SMTP_HOST_NAME = "mail.yale.edu";
# 	private final String SMTP_AUTH_USER = "";
# 	private final String SMTP_AUTH_PWD  = "";
#
# 	private final String SMTP_AUTH_PWD  = "";
# 	private final String CONTENT_TYPE = "text/html; charset=utf-8";

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

    UPLOAD_FOLDER = os.environ.get('UPLOAD_FOLDER')

    SGG_DB_CNF = os.environ.get('MYSQLDB_CNF')
    SGG_DB_HOST = os.environ.get('MYSQLDB_HOST')

    SERVER_NAME = os.environ.get('SERVER_NAME')

    PROJ_DELETE_INTERVAL = 1814400  # 3 weeks

class DevelopmentConfig(Config):
    DEBUG = True
    SECRET_KEY = os.environ.get('SECRET_KEY') or 't0p s3cr3t'
    SQLALCHEMY_DATABASE_URI = os.environ.get('DEV_DATABASE_URL') or \
        'sqlite:///' + os.path.join(basedir, 'data-dev.sqlite')
    MAIL_FLUSH_INTERVAL = 60  # one minute
    SGG_DB_NAME = os.environ.get('MYSQLDB_DB_DEV')
    SQLALCHEMY_ECHO = os.environ.get('SQLALCHEMY_ECHO') == 'True'

class TestingConfig(Config):
    TESTING = True
    SECRET_KEY = 'secret'
    SQLALCHEMY_DATABASE_URI = os.environ.get('TEST_DATABASE_URL') or \
        'sqlite:///' + os.path.join(basedir, 'data-test.sqlite')
    SGG_DB_NAME = os.environ.get('MYSQLDB_DB_TEST')


class ProductionConfig(Config):
    SQLALCHEMY_DATABASE_URI = os.environ.get('DATABASE_URL') or \
        'sqlite:///' + os.path.join(basedir, 'data.sqlite')
    SGG_DB_NAME = os.environ.get('MYSQLDB_DB')


config = {
    'development': DevelopmentConfig,
    'testing': TestingConfig,
    'production': ProductionConfig,

    'default': DevelopmentConfig
}


