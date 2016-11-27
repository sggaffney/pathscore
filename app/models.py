import os
import shutil
from datetime import datetime, timedelta
from collections import OrderedDict
from itsdangerous import TimedJSONWebSignatureSerializer as Serializer
from flask import current_app, url_for
from flask_login import current_user
from flask_security import RoleMixin, UserMixin
from . import db, login_manager
from unidecode import unidecode
from flask_security.utils import encrypt_password, verify_password

import naming_rules
from errors import ValidationError
from misc import GeneListTester, generate_random_str


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
        return generate_random_str(length=size)

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
    algorithm = db.Column(db.String(255), default='bmr_length')
    genome_size = db.Column(db.String(255), default=False)
    n_cutoff = db.Column(db.Integer)
    required_genes = db.Column(db.Text())
    ignore_genes = db.Column(db.Text())
    proj_suffix = db.Column(db.String(255))
    n_patients = db.Column(db.Integer)
    n_rejected = db.Column(db.Integer)
    n_ignored = db.Column(db.Integer)
    n_loaded = db.Column(db.Integer)
    bmr_id = db.Column(db.Integer, db.ForeignKey('bmr.bmr_id'))
    has_annot = db.Column(db.Boolean)
    has_mds = db.Column(db.Boolean)

    uploader = db.relationship('User', backref='uploads')
    bmr = db.relationship("CustomBMR", back_populates="projects")

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
            ('bmr_id', self.bmr_id),
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
                        'proj_suffix', 'bmr_id'}
        provided_keys = set(data.keys())
        unrecognized = provided_keys.difference(allowed_keys)
        if unrecognized:
            raise ValidationError("Unrecognized parameters: {}".\
                                  format(list(unrecognized)))

        algorithm = data.get('algorithm', None)
        required_genes = data.get('required_genes', None)
        ignore_genes = data.get('ignore_genes', None)
        proj_suffix = data.get('proj_suffix', None)
        bmr_id = data.get('bmr_id', None)

        if algorithm:
            if algorithm not in ['gene_count', 'gene_length', 'bmr_length']:
                raise ValidationError("Algorithm should be one of: "
                                      "gene_count, gene_length, bmr_length.")
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
        if bmr_id:
            if not bmr_id.isdigit():
                raise ValidationError("bmr_id must be a positive integer.")
            self.bmr_id = int(bmr_id)
        return self

    def load_bmr(self):
        """load bmr file into db table."""
        bmr = self.bmr  # type: CustomBMR
        if not self.bmr:
            pass  # no custom bmr, so no loading required
        bmr_path = bmr.get_path('final')
        BmrProcessor(bmr).load_final(self.file_id)
        # if not loaded:
        #     raise Exception("Failed to load bmr {!r}.".format(bmr_path))

    def remove_bmr(self):
        """delete bmr table."""
        bmr = self.bmr  # type: CustomBMR
        if not self.bmr:
            pass  # no custom bmr, so no loading required
        BmrProcessor(bmr).remove_table(self.file_id)


def create_anonymous_user():
    """Create a new anonymous user with random password."""
    temp_pswd = User.generate_random_password()
    temp_user = User(password_raw=temp_pswd, active=1)
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
    user_upload.has_annot = mut_file.has_annot
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


class CustomBMR(db.Model):
    __tablename__ = 'bmr'
    # save mut_filename, time, size, user_id in uploads table
    # with 'is_valid' field, 'run_complete'
    # run_id primary key
    # only accept this file if user has no running jobs.
    # will create folder for job: <run_id>
    bmr_id = db.Column(db.Integer, primary_key=True)
    user_id = db.Column(db.Integer, db.ForeignKey('users.id'))
    title = db.Column(db.String(255))
    tissue = db.Column(db.String(255))
    description = db.Column(db.String(255))
    upload_time = db.Column(db.DateTime(), default=datetime.utcnow)
    rand_str = db.Column(db.String(255))
    is_valid = db.Column(db.Boolean)
    n_rejected = db.Column(db.Integer)
    n_ignored = db.Column(db.Integer)
    n_loaded = db.Column(db.Integer)

    uploader = db.relationship('User', backref='bmr_files')
    projects = db.relationship("UserFile", back_populates="bmr")

    def get_filename(self, kind='final'):
        """Custom BMR filename. Kind in {'final', 'ignored', 'rejected'}."""
        assert kind in {'final', 'orig', 'ignored', 'rejected'}, "Inappropriate kind."
        keepcharacters = (' ','.', '_')
        safe_title = "".join(c for c in self.title if c.isalnum()
                              or c in keepcharacters) \
            .rstrip().replace(' ', '_')

        safe_title = unidecode(safe_title)  # returns ascii str
        # safe_suffix = ud.normalize("NFC", safe_suffix)
        filename = '_'.join([safe_title, str(self.rand_str), kind]) + '.txt'
        return filename

    def get_path(self, kind='final'):
        """Custom BMR file location. Kind in {'final', 'ignored', 'rejected'}."""
        assert kind in {'final', 'orig', 'ignored', 'rejected'}, "Inappropriate kind."
        filename = self.get_filename(kind=kind)
        dir_path = naming_rules.get_bmr_folder(self.user_id)
        return os.path.join(dir_path, filename)

    def load_final(self, proj_id):
        """Populate project-specific table with final custom bmr file."""
        file_path = self.get_path('final')
        table_name = self.get_proj_table_name(proj_id)
        cmd = u"CREATE TABLE `{table}` (`hugo_symbol` VARCHAR(255) NOT NULL, " \
              "`entrez_id` INT(11), `per_Mb` FLOAT, `length_bp` int(11), " \
              "`effective_bp` int(11) unsigned, " \
              "KEY `bmr_hugo_entrez` (`entrez_id`, `hugo_symbol`) USING HASH);"
        load_str = u"""load data local infile '{}'
            into table `{}` fields terminated by '\t'
            lines terminated by '\n' ignore 1 lines;""".format(
            file_path, table_name)
        db.session.execute(cmd.format(table=table_name))
        db.session.execute(load_str.format(table=table_name))
        db.session.commit()

    def remove_table(self, proj_id):
        """Delete final bmr table. Called after custom analyses completes."""
        table_name = self.get_proj_table_name(proj_id)
        cmd = u"drop table {};".format(table_name)
        db.session.execute(cmd)
        db.session.commit()

    def init_from_upload(self, bmr_file):
        bmr_dir = naming_rules.get_bmr_folder(self.user_id)
        if not os.path.exists(bmr_dir):
            os.mkdir(bmr_dir)
        self.rand_str = generate_random_str(6)
        bmr_file.move_file(self.get_path('orig'))
        return self

    def get_init_table_name(self):
        return 'bmr_init_{}'.format(self.user_id)

    def get_proj_table_name(self, proj_id):
        return 'bmr_{}'.format(proj_id)

    def get_url(self):
        """For API."""
        return url_for('api.get_bmr', bmr_id=self.bmr_id, _external=True)

    def export_data(self):
        """For API."""
        info_dict = OrderedDict([
            ('self_url', self.get_url()),
            ('title', self.title),
            ('tissue', self.tissue),
            ('description', self.description),
            ('upload_time', self.upload_time),
            ('n_rejected', self.n_rejected),
            ('n_unused', self.n_ignored),
            ('n_loaded', self.n_loaded)
        ])
        # url_dict = self.get_related_urls()
        # for key, val in url_dict.iteritems():
        #     info_dict[key] = val
        return info_dict

    def import_data(self, data):
        """Create CustomBMR from http form data - used by API.

        Args:
            data (dict): form data from http request.
        """
        # CHECK FOR BAD PARAMETERS
        allowed_keys = {'title', 'tissue', 'description'}
        provided_keys = set(data.keys())
        unrecognized = provided_keys.difference(allowed_keys)
        if unrecognized:
            raise ValidationError("Unrecognized parameters: {}". \
                                  format(list(unrecognized)))

        title = data.get('title', None)
        tissue = data.get('tissue', None)
        description = data.get('description', None)

        if not title:
            raise ValidationError("`title` parameter is required.")
        if len(title) > 32:
            raise ValidationError("32 character limit for `title`.")
        self.title = title
        if tissue:
            if len(tissue) > 100:
                raise ValidationError("100 char limit for `tissue`.")
            else:
                self.tissue = tissue
        if description:
            if len(tissue) > 255:
                raise ValidationError("255 char limit for `description`.")
            else:
                self.description = description

        return self


class BmrProcessor:
    """Loads data from initial file into table, given table_name and file path.

    Method initial_process performs filtering and CustomBMR object updates.
    """
    def __init__(self, bmr):
        self.bmr = bmr  # type: CustomBMR
        self.table_name = None
        self.headers = ['hugo_symbol', 'entrez_id', 'per_Mb', 'length_bp',
                        'effective_bp']

    def initial_process(self):
        """Takes CustomBMR instance; filters file; updates model."""
        loaded = False
        self.table_name = self.bmr.get_init_table_name()
        n_initial = self._populate_table()
        if n_initial:
            n_rejected = self._remove_genes_unrecognized()

        n_ignored = self._remove_genes_outwith_pathways()
        n_loaded = n_initial - n_rejected - n_ignored
        if n_loaded:
            # join to entrez_length to get lengths for scaling and extra genes
            self._fill_missing_genes()
            loaded = True
            self._save_final()
        self._update_db(loaded, n_loaded, n_rejected, n_ignored)
        db.session.execute('drop table {}'.format(self.table_name))
        db.session.commit()

    def _populate_table(self):
        """Build mutation table before filtering. Return boolean for success."""
        table_name = self.table_name
        data_path = self.bmr.get_path('orig')
        # chrom?: `chrom` VARCHAR(255) DEFAULT NULL,
        # `effective_bp` FLOAT NOT NULL,
        create_str = u"""CREATE TABLE `{}` (
          `hugo_symbol` VARCHAR(255) NOT NULL,
          `entrez_id` INT(11) DEFAULT NULL,
          `per_Mb` FLOAT DEFAULT NULL,
          KEY `bmr_hugo_entrez` (`entrez_id`, `hugo_symbol`) USING HASH);"""\
            .format(table_name)
        load_str = u"""load data local infile '{}'
            into table `{}` fields terminated by '\t'
            lines terminated by '\n' ignore 1 lines;""".format(
            data_path, table_name)
        cmd = u"select count(*) from `{}` m;".format(self.table_name)
        db.session.execute(create_str)
        db.session.execute(load_str)
        n_initial = db.session.execute(cmd).scalar()
        db.session.commit()
        return n_initial

    def _remove_genes_unrecognized(self):
        cmd1 = u"""SELECT m.* FROM `{}` m
            LEFT JOIN refs.ncbi_entrez n ON m.entrez_id = n.geneId
            WHERE n.geneId IS NULL OR m.hugo_symbol <> n.symbol;""".format(
            self.table_name)
        cmd2 = u"""delete from m using `{}` m
          LEFT JOIN refs.ncbi_entrez n ON m.entrez_id = n.geneId
            WHERE n.geneId IS NULL OR m.hugo_symbol <> n.symbol;""" \
            .format(self.table_name)
        # EXPORT REJECTED GENES
        result = db.session.execute(cmd1)
        n_rejected = result.rowcount
        if n_rejected > 0:
            self._save_mutations_subset(result, self.bmr.get_path(kind='rejected'))
            # DELETE EXTRA GENE LINES
            db.session.execute(cmd2)
            db.session.commit()
        return n_rejected

    def _remove_genes_outwith_pathways(self, rejected_path=None):
        """
        Remove genes outwith pathways, save to reject file,
        return remaining mutation count.

        return: remaining mutation count
        :rtype : int
        """
        cmd1 = """SELECT m.* FROM `{}` m
          LEFT JOIN (SELECT DISTINCT entrez_id FROM refs.pathway_gene_link) l
          ON m.entrez_id = l.entrez_id WHERE l.entrez_id IS NULL;""" \
            .format(self.table_name)
        cmd2 = u"""delete from m using `{}` m
          LEFT JOIN (SELECT DISTINCT entrez_id FROM refs.pathway_gene_link) l
          ON m.entrez_id = l.entrez_id WHERE l.entrez_id IS NULL;""" \
            .format(self.table_name)
        # EXPORT EXTRA GENE LINES
        result = db.session.execute(cmd1)
        n_ignored = result.rowcount
        if n_ignored > 0:
            self._save_mutations_subset(result,
                                        self.bmr.get_path(kind='ignored'))
            # DELETE EXTRA GENE LINES
            db.session.execute(cmd2)
            db.session.commit()
        return n_ignored

    @staticmethod
    def _save_mutations_subset(result, save_path):
        """Export mutation subset."""
        with open(save_path, 'w') as out:
            for row in result:
                out.write('\t'.join([str(i) for i in row]) + '\n')

    def _fill_missing_genes(self):
        """Combine with refs.entrez_length for length_bp; get effective_bp."""
        cmd_alter = u"ALTER TABLE {table} ADD COLUMN length_bp INT(11), " \
                    "ADD COLUMN effective_bp INT(11);"
        cmd_fetch_len = u"UPDATE refs.entrez_length e LEFT JOIN {table} b " \
                        "ON e.`entrez_id` = b.`entrez_id` " \
                        "SET b.length_bp = e.length_bp;"
        cmd_extra = u"INSERT INTO {table} (hugo_symbol, entrez_id, per_Mb, " \
                    "length_bp, effective_bp) SELECT e.hugo_symbol, " \
                    "e.entrez_id, NULL, e.length_bp, NULL FROM " \
                    "refs.entrez_length e LEFT JOIN {table} b " \
                    "ON e.`entrez_id` = b.`entrez_id` " \
                    "WHERE b.`entrez_id` IS NULL;"
        cmd_permb = u"UPDATE {table} b INNER JOIN (SELECT AVG(per_Mb) AS " \
                    "per_Mb FROM {table}) a SET b.per_Mb = a.per_Mb " \
                    "WHERE b.per_Mb IS NULL;"
        cmd_effective = u"UPDATE {table} SET `effective_bp` = `per_Mb` * " \
                        u"`length_bp`;"
        db.session.execute(cmd_alter.format(table=self.table_name))
        db.session.execute(cmd_fetch_len.format(table=self.table_name))
        db.session.execute(cmd_extra.format(table=self.table_name))
        db.session.execute(cmd_permb.format(table=self.table_name))
        db.session.execute(cmd_effective.format(table=self.table_name))
        db.session.commit()

    def _update_db(self, loaded, n_loaded, n_rejected, n_ignored):
        """create db entry, even if there are no valid/kept genes"""
        self.bmr.is_valid = loaded
        self.bmr.n_rejected = n_rejected
        self.bmr.n_ignored = n_ignored
        self.bmr.n_loaded = n_loaded
        db.session.add(self.bmr)
        db.session.commit()

    def _save_final(self):
        """Save filtered file to final path."""
        # bmr knows final path.
        tmp_dir = current_app.config['TEMP_FOLDER']
        user_id = self.bmr.user_id
        time = datetime.utcnow().strftime('%Y-%m-%d_%H%M%S_%f')
        rand = generate_random_str(3)
        tmp_name = 'bmr_{user}_{time}_{rand}.txt'.format(
            user=user_id, time=time, rand=rand)
        tmp_path = os.path.join(tmp_dir, tmp_name)
        final_path = self.bmr.get_path('final')
        columns_str = ', '.join([repr(i) for i in self.headers])
        cmd_fetch = u"select {columns} union all select * from {table} " \
                    u"into outfile {tmp!r};".\
            format(columns=columns_str, table=self.table_name, tmp=tmp_path)
        db.session.execute(cmd_fetch)
        shutil.move(tmp_path, final_path)
