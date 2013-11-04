
from distutils.core import setup
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
ext_modules = [Extension("admixprop", ["admixprop.pyx", "admixprop.c"]), \
       Extension("allelefreq", ["allelefreq.pyx", "allelefreq.c"], \
        library_dirs=["/usr/local/lib/"], \
        libraries=["gsl","gslcblas"], \
        extra_compile_args=["-g"], extra_link_args=["-g"]),
       Extension("marglikehood", ["marglikehood.pyx", "marglikehood.c"])]

setup(
    name = 'variables',
    cmdclass = {'build_ext': build_ext},
    include_dirs=[numpy.get_include(), '.'],
    ext_modules = ext_modules
)
