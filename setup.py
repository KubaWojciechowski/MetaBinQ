from distutils.core import setup, Extension

module1 = Extension('tnCounter',
                    sources = ['tnCounter.c'])

setup (name = 'tnCounter',
       version = '1.0',
       description = 'This package counts tetra nucelotides',
       ext_modules = [module1])
