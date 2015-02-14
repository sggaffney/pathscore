import unicodedata as ud
from subprocess import check_call
from glob import glob
import os


def simplify_string(odd_string):
    keepcharacters = (' ', '.', '_')
    safe_string = "".join(c for c in odd_string if c.isalnum()
                          or c in keepcharacters)\
        .rstrip().replace(' ', '_')
    if safe_string:
        safe_string = safe_string.encode('latin-1', 'ignore')\
            .decode('latin-1')
        safe_string = ud.normalize("NFC", safe_string)
    return safe_string


def zip_svgs(proj_dir):
    """Zip and delete all svgs in matrix_svg, pathways_svg."""
    svg_dirs = ['pathways_svg', 'matrix_svg']
    for plot_dir in svg_dirs:
        file_names = glob(os.path.join("{}".format(proj_dir), plot_dir, '*.svg'))
        for file_name in file_names:
            check_call(['/usr/local/bin/svgo', file_name])
            check_call(['gzip', file_name])
            os.rename(file_name + '.gz', file_name + 'z')