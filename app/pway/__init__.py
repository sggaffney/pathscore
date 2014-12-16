from flask import Blueprint

pway = Blueprint('pway', __name__)


class FileTester():
    want_headers = ['hugo_symbol', 'entrez_id', 'patient_id',
                    'variant_classification']

    def __init__(self, file_path):
        self.good_headers = False  # will overwrite this if headers are valid
        self.data_present = False  # will overwrite this if line after headers

        self._check_headers_data(file_path)

    def _check_headers_data(self, file_path):
        line = None
        with open(file_path, 'rU') as file:
            for line in file:
                if line.startswith('#') or line.strip() == '':
                    continue
                else:
                    break
            # want at least 2 lines now. 1 for header, 1 for data.
            if line:
                headers = [i.lower() for i in line.strip('\n').split('\t')]
                if FileTester._compare_headers(headers):
                    self.good_headers = True
            for line in file:
                if line.strip('\n'):  # data line present
                    self.data_present = True

    @staticmethod
    def _compare_headers(headers):
        """Return True if provided headers list matches desired headers."""
        if len(FileTester.want_headers) != len(headers):
            return False
        for pair in zip(FileTester.want_headers, headers):
            if pair[0] != pair[1]:
                return False
        return True

from . import routes


