from distutils.core import setup, Extension
import os

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


module1 = Extension('clapy.ccla',
                    define_macros = [('MAJOR_VERSION', '1'),
                                     ('MINOR_VERSION', '0')],
                    language = "c++",
                    extra_compile_args=['-std=c++11','-I/usr/include/gsl'],
                    extra_link_args=['-lgsl','-lgslcblas'],
                    sources = ['src/ccla.cpp',
                        'src/types.cpp',
                        'src/cell.cpp',
                        'src/rand.cpp'])

setup (name = 'clapy',
       version = '1.0',
       description = 'This is a demo package',
       packages = ['clapy'],
       author = '',
       author_email = '',
       url = '',
       long_description=read('README'),
       ext_modules = [module1],
       )
