import os
from datetime import datetime, timedelta
import random
import string
from collections import OrderedDict
from itsdangerous import TimedJSONWebSignatureSerializer as Serializer
from flask import current_app, url_for
from flask_login import current_user
from flask_security import RoleMixin, UserMixin
from . import db, login_manager
from unidecode import unidecode
from flask.ext.security.utils import encrypt_password, verify_password
import naming_rules
from errors import ValidationError
from misc import GeneListTester

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
                      nullable=True, unique=True, index=True)
    password = db.Column(db.String(255))
    member_since = db.Column(db.DateTime(), default=datetime.utcnow)
    active = db.Column(db.Boolean())
    confirmed_at = db.Column(db.DateTime())
    roles = db.relationship('Role', secondary=roles_users,
                            backref=db.backref('users', lazy='dynamic'))

    @property
    def password_raw(self):
        raise AttributeError('password is not a readable attribute')

    @password_raw.setter
    def password_raw(self, password):
        """Used by manager adduser to has password."""
        self.password = encrypt_password(password)

    def verify_password(self, password):
        return verify_password(password, self.password)
        # return check_password_hash(self.password_hash, password)

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

    @staticmethod
    def generate_random_password(size=12):
        # via http://stackoverflow.com/a/2257449
        chars = string.ascii_uppercase + string.digits + string.ascii_lowercase
        return ''.join(random.SystemRandom().choice(chars) for _ in range(size))

    @staticmethod
    def get_guest_username(int_id):
        return 'guest_' + hex(int(int_id))[2:]

    def get_expiry_date(self):
        """If user is anonymous, return deletion date, else return None."""
        if 'anonymous' in [r.name for r in self.roles]:
            max_age_days = current_app.config['ANONYMOUS_MAX_AGE_DAYS']
            return self.member_since + timedelta(days=max_age_days)
        else:
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
    # with 'is_valid' field, 'run_complete'
    # run_id primary key
    # only accept this file if user has no running jobs.
    # will create folder for job: <run_id>
    file_id = db.Column(db.Integer, primary_key=True)
    user_id = db.Column(db.Integer, db.ForeignKey('users.id'))
    filename = db.Column(db.String(255))
    upload_time = db.Column(db.DateTime(), default=datetime.utcnow)
    is_valid = db.Column(db.Boolean)
    is_queued = db.Column(db.Boolean)
    run_complete = db.Column(db.Boolean)
    algorithm = db.Column(db.String(255), default='gene_length')
    genome_size = db.Column(db.String(255), default=False)
    n_cutoff = db.Column(db.Integer)
    required_genes = db.Column(db.Text())
    ignore_genes = db.Column(db.Text())
    proj_suffix = db.Column(db.String(255))
    uploader = db.relationship('User', backref='uploads')
    n_patients = db.Column(db.Integer)
    n_rejected = db.Column(db.Integer)
    n_ignored = db.Column(db.Integer)
    n_loaded = db.Column(db.Integer)

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

    @property
    def ignore_short(self):
        ignore_list = self.ignore_genes.split(',')
        ignore_str = ','.join(ignore_list[0:5])
        if len(ignore_list) > 5:
            ignore_str += '...'
        return ignore_str

    def get_expiry_date(self):
        """Return project expiry date.

        This is user expiry date for anonymous users, standard expiry date for
        general users, and None for vips.
        ."""
        user_expiry = self.uploader.get_expiry_date()
        if user_expiry:
            return user_expiry
        if 'vip' in [r.name for r in self.roles]:
            return None
        timediff = timedelta(days=current_app.config['PROJ_MAX_AGE_DAYS'])
        return self.upload_time + timediff

    def get_status(self):
        """Status message for views."""
        if self.is_queued:
            return "Queued"
        if self.run_complete:
            return "Complete"
        if self.run_complete is None:
            return "Failed"
        if not self.run_complete:
            return "Running"
        return "Unknown"

    def get_url(self):
        return url_for('api.get_project', file_id=self.file_id, _external=True)

    def get_related_urls(self):
        url_dict = OrderedDict()
        if self.run_complete:
            url_dict['status'] = 'Project complete.'
            url_dict['flat_url'] = url_for('pway.results', proj=self.file_id,
                                           _external=True)
            url_dict['scatter_url'] = url_for('pway.scatter', proj=self.file_id,
                                              _external=True)
            url_dict['tree_url'] = url_for('pway.tree', proj=self.file_id,
                                           _external=True)
        elif self.run_complete == 0:
            url_dict['status'] = 'Not yet complete.'
        elif self.run_complete is None:
            url_dict['status'] = 'Failed.'
        url_dict['archive_url'] = url_for('api.archive', proj=self.file_id,
                                          _external=True)
        url_dict['filtered_unused'] = url_for('pway.get_filtered',
                                              proj=self.file_id,
                                              type='ignored', _external=True)
        url_dict['filtered_rejected'] = url_for('pway.get_filtered',
                                                proj=self.file_id,
                                                type='rejected', _external=True)
        return url_dict

    def export_data(self):
        info_dict = OrderedDict([
            ('self_url', self.get_url()),
            ('name', self.get_local_filename()),
            ('upload_time', self.upload_time),
            ('algorithm', self.algorithm),
            ('required_genes', self.required_genes),
            ('ignore_genes', self.ignore_genes),
            ('n_patients', self.n_patients),
            ('n_rejected', self.n_rejected),
            ('n_unused', self.n_ignored),
            ('n_loaded', self.n_loaded)
        ])
        url_dict = self.get_related_urls()
        for key, val in url_dict.iteritems():
            info_dict[key] = val
        return info_dict

    def import_data(self, data):
        """Create UserFile from http form data - used by API.

        Args:
            data (dict): form data from http request.
        """
        # CHECK FOR BAD PARAMETERS
        allowed_keys = {'algorithm', 'required_genes', 'ignore_genes',
                        'proj_suffix'}
        provided_keys = set(data.keys())
        unrecognized = provided_keys.difference(allowed_keys)
        if unrecognized:
            raise ValidationError("Unrecognized parameters: {}".\
                                  format(list(unrecognized)))

        algorithm = data.get('algorithm', None)
        required_genes = data.get('required_genes', None)
        ignore_genes = data.get('ignore_genes', None)
        proj_suffix = data.get('proj_suffix', None)

        if algorithm:
            if algorithm not in ['gene_count', 'gene_length']:
                raise ValidationError("Algorithm should be one of: "
                                      "gene_count, gene_length.")
            else:
                self.algorithm = algorithm
        if required_genes:
            if not GeneListTester.is_valid(required_genes):
                raise ValidationError("Required genes must comma-separated, "
                                      "no spaces.")
            elif len(required_genes) > 600:
                raise ValidationError("600 char limit for required_genes.")
            else:
                self.required_genes = required_genes
        if ignore_genes:
            if not GeneListTester.is_valid(ignore_genes):
                raise ValidationError("Ignored genes must comma-separated, "
                                      "no spaces.")
            elif len(ignore_genes) > 600:
                raise ValidationError("600 char limit for ignore_genes.")
            else:
                self.ignore_genes = ignore_genes
        if proj_suffix:
            if len(proj_suffix) > 40:
                raise ValidationError("40 char limit for proj_suffix.")
            self.proj_suffix = proj_suffix
        return self


def create_anonymous_user():
    """Create a new anonymous user with random password."""
    temp_pswd = User.generate_random_password()
    temp_user = User(password_raw=temp_pswd)
    db.session.add(temp_user)
    db.session.commit()
    temp_uname = User.get_guest_username(temp_user.id)
    temp_user.email = temp_uname
    temp_user.confirmed_at = datetime.utcnow()
    temp_user.roles.append(Role.query.filter_by(name='anonymous').one())
    db.session.commit()
    return temp_user, temp_pswd


def initialize_project(user_upload=None, mut_file=None):
    """Tweak upload object. Initialize folder and move mutation file.

    Args:
        user_upload (UserFile): upload Model instance
        mut_file (MutationFile): mutation file info and methods
    """

    # CREATE USERFILE FROM FORM AND ADD TO DB
    user_upload.is_valid = True
    user_upload.run_complete = False
    user_upload.is_queued = True
    db.session.add(user_upload)
    db.session.commit()

    current_app.logger.info("Project {} uploaded by {}".format(
        user_upload.file_id, current_user.email))

    # MAKE PROJECT FOLDER
    user_folder = naming_rules.get_user_folder(user_upload.user_id)
    proj_folder = naming_rules.get_project_folder(user_upload)
    if not os.path.exists(user_folder):
        os.mkdir(user_folder)
    os.mkdir(proj_folder)

    # MOVE TEMP FILE
    file_path = os.path.join(proj_folder, user_upload.get_local_filename())
    mut_file.move_file(file_path)  # move to project folder

    return user_upload, proj_folder, file_path
