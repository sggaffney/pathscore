from datetime import datetime
from werkzeug.security import generate_password_hash, check_password_hash
from itsdangerous import TimedJSONWebSignatureSerializer as Serializer
from flask import request, current_app
from flask_login import UserMixin
from flask_security import RoleMixin
from . import db, login_manager
import unicodedata as ud
from unidecode import unidecode
from flask_wtf.file import FileField, FileAllowed, FileRequired
from flask.ext.security.utils import encrypt_password

# Define models
roles_users = db.Table('roles_users',
                       db.Column('user_id', db.Integer(),
                                 db.ForeignKey('users.id')),
                       db.Column('role_id', db.Integer(),
                                 db.ForeignKey('roles.id')))


class Role(db.Model, RoleMixin):
    __tablename__ = 'roles'
    id = db.Column(db.Integer(), primary_key=True)
    name = db.Column(db.String(80), unique=True)
    description = db.Column(db.String(255))
    uploads_pw = db.Column(db.Integer())


class User(UserMixin, db.Model):
    __tablename__ = 'users'
    id = db.Column(db.Integer, primary_key=True)
    email = db.Column(db.String(64),
                      nullable=False, unique=True, index=True)
    is_admin = db.Column(db.Boolean)
    password = db.Column(db.String(255))
    name = db.Column(db.String(64))
    member_since = db.Column(db.DateTime(), default=datetime.utcnow)
    confirmed_at = db.Column(db.DateTime())
    active = db.Column(db.Boolean())
    roles = db.relationship('Role', secondary=roles_users,
                            backref=db.backref('users', lazy='dynamic'))

    @property
    def password_raw(self):
        raise AttributeError('password is not a readable attribute')

    @password_raw.setter
    def password_raw(self, password):
        """Used by manager adduser to has password."""
        self.password = encrypt_password(password)

    # def verify_password(self, password):
    #     return check_password_hash(self.password_hash, password)

    def get_api_token(self, expiration=300):
        s = Serializer(current_app.config['SECRET_KEY'], expiration)
        return s.dumps({'user': self.id}).decode('utf-8')

    @staticmethod
    def validate_api_token(token):
        s = Serializer(current_app.config['SECRET_KEY'])
        try:
            data = s.loads(token)
        except:
            return None
        id = data.get('user')
        if id:
            return User.query.get(id)
        return None

    def get_auth_token(self, expiration=300):
        s = Serializer(current_app.config['SECRET_KEY'], expiration)
        return s.dumps({'user': self.id}).decode('utf-8')

    @staticmethod
    def validate_auth_token(token):
        s = Serializer(current_app.config['SECRET_KEY'])
        try:
            data = s.loads(token)
        except:
            return None
        id = data.get('user')
        if id:
            return User.query.get(id)
        return None

@login_manager.user_loader
def load_user(user_id):
    return User.query.get(int(user_id))

@login_manager.token_loader
def load_token(token):
    s = Serializer(current_app.config['SECRET_KEY'])
    try:
        data = s.loads(token)
    except:
        return None
    id = data.get('user')
    if id:
        return User.query.get(id)
    return None


class UserFile(db.Model):
    __tablename__ = 'uploads'
    # save mut_filename, time, size, user_id in uploads table
    # with 'is_valid' field, 'run_accepted', 'run_complete'
    # run_id primary key
    # only accept this file if user has no running jobs.
    # will create folder for job: <run_id>
    file_id = db.Column(db.Integer, primary_key=True)
    user_id = db.Column(db.Integer, db.ForeignKey('users.id'))
    filename = db.Column(db.String(255))
    upload_time = db.Column(db.DateTime(), default=datetime.utcnow)
    is_valid = db.Column(db.Boolean)
    run_accepted = db.Column(db.Boolean)
    run_complete = db.Column(db.Boolean)
    genome_size = db.Column(db.String(255), default=False)
    n_cutoff = db.Column(db.Integer)
    required_genes = db.Column(db.Text())
    ignore_genes = db.Column(db.Text())
    proj_suffix = db.Column(db.String(255))
    uploader = db.relationship('User', backref='uploads')

    def get_local_filename(self):
        keepcharacters = (' ', '.', '_')
        safe_suffix = None
        if self.proj_suffix:
            safe_suffix = "".join(c for c in self.proj_suffix if c.isalnum()
                                  or c in keepcharacters)\
                .rstrip().replace(' ', '_')
        if safe_suffix:
            safe_suffix = unidecode(safe_suffix)  # returns ascii str
            # safe_suffix = ud.normalize("NFC", safe_suffix)
            use_title = safe_suffix  # str(safe_suffix)
        else:
            use_title = self.filename
        return '_'.join([str(self.file_id), use_title])

    def get_table_name(self):
        temp_id = str(self.file_id)
        return 'mutations_{}'.format(temp_id)

