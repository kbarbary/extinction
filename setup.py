#!/usr/bin/env python
import os
from setuptools import setup
from setuptools.extension import Extension
import re

import numpy
from Cython.Build import cythonize

fname = "extinction.pyx"

sourcefiles = [fname, os.path.join("extern", "bs.c")]
dependsfiles = [os.path.join("extern", "bs.h"), os.path.join("extern", "bsplines.pxi")]
include_dirs = [numpy.get_include(), "extern"]
extensions = [
    Extension(
        "extinction",
        sourcefiles,
        include_dirs=include_dirs,
        depends=dependsfiles,
        extra_compile_args=["-std=c99"],
    )
]

extensions = cythonize(extensions)

# Synchronize version from code.
version = re.findall(r"__version__ = \"(.*?)\"", open(fname).read())[0]

classifiers = [
    "Development Status :: 4 - Beta",
    "Programming Language :: Python :: 3",
    "License :: OSI Approved :: MIT License",
    "Topic :: Scientific/Engineering :: Astronomy",
    "Intended Audience :: Science/Research",
]

setup(
    name="extinction",
    version=version,
    description="Fast interstellar dust extinction laws",
    long_description="documentation: http://extinction.readthedocs.io",
    license="MIT",
    classifiers=classifiers,
    url="http://github.com/kbarbary/extinction",
    author="Kyle Barbary",
    author_email="kylebarbary@gmail.com",
    ext_modules=extensions,
    install_requires=["numpy>=1.13.3"],
)
