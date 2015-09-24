import os
import hashlib
from datetime import datetime
from flask import Blueprint, current_app, g
from app.get_effective_pathways import path_info_dict

pway = Blueprint('pway', __name__)

class TempFile:
    def __init__(self, filestorage):
        self.filestorage = filestorage
        timestr = datetime.utcnow().strftime('%Y-%m-%d_%H%M%S_%f_')
        md5str = self._generate_stream_md5()
        self.temp_name = timestr + md5str
        self.path = self._create_temp_copy()

    def _generate_stream_md5(self, blocksize=2**20):
        m = hashlib.md5()
        for chunk in iter(lambda: self.filestorage.stream.read(blocksize), b''):
            m.update( chunk )
        return m.hexdigest()

    def _create_temp_copy(self):
        temp_dir = current_app.config['TEMP_FOLDER']
        temp_path = os.path.join(temp_dir, self.temp_name)
        self.filestorage.seek(0)
        self.filestorage.save(temp_path)
        return temp_path


class FileTester:
    want_headers = ['hugo_symbol', 'entrez_id', 'patient_id']
    int_columns = [1]

    def __init__(self, file_path):
        self.data_issues = []
        self.line_endings = None
        self._check_headers_data(file_path)  # get data_issues and line_endings

    def _check_headers_data(self, file_path):
        ind = -1
        with open(file_path, 'rU') as file:
            line = file.readline()
            # want at least 2 lines now. 1 for header, 1 for data.
            headers = [i.lower() for i in line.strip('\n').split('\t')]
            if not FileTester._compare_headers(headers):
                self.data_issues.append("Invalid headers.")
                return
            if not line:
                self.data_issues.append("File is empty.")
                return
            for ind, line in enumerate(file):
                if not self._line_valid(line):
                    self.data_issues.append("Line {} is invalid".format(
                        ind + 1))
                    return
            if not ind + 1:
                self.data_issues.append("No data found.")
                return
            self.line_endings = file.newlines
            # could add minimum line count here
            # if ind < 9:
            #     self.data_issues.append("Need more data.")
            #     return

    @staticmethod
    def _line_valid(line):
        vals = line.strip('\n').split('\t')
        if len(vals) != len(FileTester.want_headers):
            return False
        int_vals = [vals[j] for j in FileTester.int_columns]
        if False in [i.isdigit() for i in int_vals]:
            print line
            return False
        return True

    @staticmethod
    def _compare_headers(headers):
        """Return True if provided headers list matches desired headers."""
        if len(FileTester.want_headers) != len(headers):
            return False
        for pair in zip(FileTester.want_headers, headers):
            if pair[0] != pair[1]:
                return False
        return True

@pway.context_processor
def inject_n_pathways():
    return dict(n_pathways=len(path_info_dict))

@pway.before_request
def get_n_pathways():
    g.n_pathways = len(path_info_dict)

from . import routes


