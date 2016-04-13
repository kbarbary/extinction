#!/usr/bin/env python
import os
from distutils.core import setup
from distutils.extension import Extension
import re

if os.path.exists("extinction.pyx"):
    USE_CYTHON = True
    fname = "extinction.pyx"
else:
    USE_CYTHON = False
    fname = "extinction.c"

extensions = [Extension("extinction", [fname])]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

# Synchronize version from code.
version = re.findall(r"__version__ = \"(.*?)\"", open(fname).read())[0]

classifiers = []

setup(name="extinction", 
      version=version,
      description="",
      long_description="",
      license="MIT",
      classifiers=classifiers,
      url="",
      author="Kyle Barbary",
      author_email="kylebarbary@gmail.com",
      ext_modules=extensions,
      requires=["numpy"])
