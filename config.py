import os
basedir = os.path.abspath(os.path.dirname(__file__))


# 	private final String SMTP_HOST_NAME = "mail.yale.edu";
# 	private final String SMTP_AUTH_USER = "";
# 	private final String SMTP_AUTH_PWD  = "";
#
# 	private final String SMTP_AUTH_PWD  = "";
# 	private final String CONTENT_TYPE = "text/html; charset=utf-8";


class Config:
    SECRET_KEY = os.environ.get('SECRET_KEY')

    MAIL_SERVER = os.environ.get('MAIL_SERVER')
    MAIL_PORT = int(os.environ.get('MAIL_PORT'))
    MAIL_USE_TLS = bool(os.environ.get('MAIL_USE_TLS') == 'True') or False
    MAIL_USE_SSL = bool(os.environ.get('MAIL_USE_SSL') == 'True') or False
    MAIL_USERNAME = os.environ.get('MAIL_USERNAME')
    MAIL_PASSWORD = os.environ.get('MAIL_PASSWORD')
    MAIL_DEFAULT_SENDER = os.environ.get('MAIL_DEFAULT_SENDER') or ''
    MAIL_SENDER = os.environ.get('MAIL_SENDER')
    # MAIL_FLUSH_INTERVAL = 3600  # one hour
    MAIL_ERROR_RECIPIENT = os.environ.get('MAIL_ERROR_RECIPIENT')

    SECURITY_EMAIL_SENDER = os.environ.get('SECURITY_EMAIL_SENDER')
    SECURITY_PASSWORD_HASH = os.environ.get('SECURITY_PASSWORD_HASH')
    SECURITY_PASSWORD_SALT = os.environ.get('SECURITY_PASSWORD_SALT')

    UPLOAD_FOLDER = os.environ.get('UPLOAD_FOLDER')

    SGG_DB_CNF = os.environ.get('MYSQLDB_CNF')
    SGG_DB_HOST = os.environ.get('MYSQLDB_HOST')


class DevelopmentConfig(Config):
    DEBUG = True
    SECRET_KEY = os.environ.get('SECRET_KEY') or 't0p s3cr3t'
    SQLALCHEMY_DATABASE_URI = os.environ.get('DEV_DATABASE_URL') or \
        'sqlite:///' + os.path.join(basedir, 'data-dev.sqlite')
    MAIL_FLUSH_INTERVAL = 60  # one minute
    SGG_DB_NAME = os.environ.get('MYSQLDB_DB_DEV')


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


