
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy
import sys

# setup bed parser
ext_modules = [Extension("parse_bed", ["parse_bed.pyx"])]

setup(
    name = 'parse_bed',
    cmdclass = {'build_ext': build_ext},
    include_dirs=[numpy.get_include(), '.'],
    ext_modules = ext_modules
)

# setup structure format parser
ext_modules = [Extension("parse_str", ["parse_str.pyx"])]

setup(
    name = 'parse_str',
    cmdclass = {'build_ext': build_ext},
    include_dirs=[numpy.get_include(), '.'],
    ext_modules = ext_modules
)

# setup fastStructure
ext_modules = [Extension("fastStructure", ["fastStructure.pyx"])]

setup(
    name = 'fastStructure',
    cmdclass = {'build_ext': build_ext},
    include_dirs=[numpy.get_include(), '.', 'vars/'],
    ext_modules = ext_modules
)
