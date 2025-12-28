# Setup script for building Cython extensions
# Run with: python setup_cython.py build_ext --inplace

from setuptools import setup
from Cython.Build import cythonize
import numpy as np

setup(
    ext_modules=cythonize(
        "app/comb_functions.pyx",
        compiler_directives={'language_level': "3"}
    ),
    include_dirs=[np.get_include()],
)
