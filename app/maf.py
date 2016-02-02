import os
from datetime import datetime
import hashlib
from flask import current_app
from errors import ValidationError


class MutationFile(object):
    """Holds file data, performs file tests, saves to appropriate locations."""

    want_headers = ['hugo_symbol', 'entrez_id', 'patient_id']
    int_columns = [1]

    def __init__(self, filestorage):
        """Save data to temporary file."""
        self.filestorage = filestorage
        self.temp_path = self._save_temp_file()
        self.line_endings = None
        # self.file_path = file_path
        self._check_headers_data()  # check for issues, get line_endings

    def _save_temp_file(self):
        timestr = datetime.utcnow().strftime('%Y-%m-%d_%H%M%S_%f_')
        md5str = self._generate_stream_md5()
        temp_name = timestr + md5str
        temp_dir = current_app.config['TEMP_FOLDER']
        temp_path = os.path.join(temp_dir, temp_name)
        self.filestorage.seek(0)
        self.filestorage.save(temp_path)
        return temp_path

    def _generate_stream_md5(self, blocksize=2**20):
        m = hashlib.md5()
        for chunk in iter(lambda: self.filestorage.stream.read(blocksize), b''):
            m.update(chunk)
        return m.hexdigest()

    def _check_headers_data(self):
        ind = -1
        with open(self.temp_path, 'rU') as tempfile:
            line = tempfile.readline()
            # want at least 2 lines now. 1 for header, 1 for data.
            headers = [i.lower() for i in line.strip('\n').split('\t')]
            if not self._compare_headers(headers):
                raise ValidationError("Invalid mutation file headers.")
            if not line:
                raise ValidationError("Mutation file is empty.")
            for ind, line in enumerate(tempfile):
                if not self._line_valid(line):
                    raise ValidationError("Line {} is invalid".format(ind + 1))
            if not ind + 1:
                raise ValidationError("No data found in mutation file.")
            self.line_endings = tempfile.newlines
            # could add minimum line count here
            # if ind < 9:
            #     raise ValidationError("Need more data.")

    def _line_valid(self, line):
        vals = line.strip('\n').split('\t')
        if len(vals) != len(self.want_headers):
            return False
        int_vals = [vals[j] for j in self.int_columns]
        if False in [i.isdigit() for i in int_vals]:
            print line
            return False
        return True

    def _compare_headers(self, headers):
        """Return True if provided headers list matches desired headers."""
        if len(self.want_headers) != len(headers):
            return False
        for pair in zip(self.want_headers, headers):
            if pair[0] != pair[1]:
                return False
        return True

    def move_file(self, file_path):
        """Move validated file to new path."""
        if self.line_endings == '\n':
            os.rename(self.temp_path, file_path)
        else:  # rewrite file with \n line endings
            with open(self.temp_path, 'rU') as tempfile:
                with open(file_path, 'w') as out:
                    for line in tempfile:
                        out.write(line)
            os.remove(self.temp_path)
