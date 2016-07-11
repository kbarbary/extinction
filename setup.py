#!/usr/bin/env python
import os
from setuptools import setup
from setuptools.extension import Extension
import re

import numpy

if os.path.exists("extinction.pyx"):
    USE_CYTHON = True
    fname = "extinction.pyx"
else:
    USE_CYTHON = False
    fname = "extinction.c"

include_dirs = [numpy.get_include()]
extensions = [Extension("extinction", [fname], include_dirs=include_dirs)]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

# Synchronize version from code.
version = re.findall(r"__version__ = \"(.*?)\"", open(fname).read())[0]

classifiers = [
    "Development Status :: 3 - Alpha",
    "Programming Language :: Python :: 2",
    "Programming Language :: Python :: 3",
    "License :: OSI Approved :: MIT License",
    "Topic :: Scientific/Engineering :: Astronomy",
    "Intended Audience :: Science/Research"]

setup(name="extinction", 
      version=version,
      description="Fast interstellar dust extinction laws",
      long_description="documentation: http://extinction.readthedocs.io",
      license="MIT",
      classifiers=classifiers,
      url="http://github.com/kbarbary/extinction",
      author="Kyle Barbary",
      author_email="kylebarbary@gmail.com",
      ext_modules=extensions,
      install_requires=["numpy", "scipy"])
