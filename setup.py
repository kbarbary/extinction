#!/usr/bin/env python
import os
import re
import sys

import numpy
from setuptools import setup
from setuptools.extension import Extension


classifiers = [
    "Development Status :: 4 - Beta",
    "Programming Language :: Python :: 3",
    "License :: OSI Approved :: MIT License",
    "Topic :: Scientific/Engineering :: Astronomy",
    "Intended Audience :: Science/Research",
]

# Synchronize version from code.
fname = "extinction.pyx"
version = re.findall(r"__version__ = \"(.*?)\"", open(fname).read())[0]

# Build Cython extension
source_files = [fname, os.path.join("extern", "bs.c")]
depends_files = [
    os.path.join("extern", "bs.h"),
    os.path.join("extern", "bsplines.pxi")
]
include_dirs = [numpy.get_include(), "extern"]
extensions = [
    Extension(
        "extinction",
        source_files,
        include_dirs=include_dirs,
        depends=depends_files,
        extra_compile_args=["-std=c99"],
    )
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
    python_requires=">=3.5",
)
