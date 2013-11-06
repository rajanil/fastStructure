
from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy
import sys

# setup utility functions
ext_modules = [Extension("utils", ["utils.pyx"])]

setup(
    name = 'utility functions',
    cmdclass = {'build_ext': build_ext},
    include_dirs=[numpy.get_include(), '.'],
    ext_modules = ext_modules
)

# setup variable updates 
ext_modules = [Extension("admixprop", ["admixprop.pyx", "C_admixprop.c"]), \
               Extension("allelefreq", ["allelefreq.pyx", "C_allelefreq.c"], \
                libraries=["gsl","gslcblas"]), \
               Extension("marglikehood", ["marglikehood.pyx", "C_marglikehood.c"])]
ext_modules = cythonize(ext_modules)

setup(
    name = 'variables',
    author = 'Anil Raj',
    version = '1.0',
    author_email = 'rajanil@stanford.edu',
    cmdclass = {'build_ext': build_ext},
    include_dirs=[numpy.get_include(), '.'],
    ext_modules = ext_modules
)
