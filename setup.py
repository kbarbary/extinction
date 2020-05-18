#!/usr/bin/env python
import os
import re
import sys
from setuptools import setup
from setuptools.extension import Extension
from setuptools.dist import Distribution


if os.path.exists("extinction.pyx"):
    USE_CYTHON = True
    fname = "extinction.pyx"
else:
    USE_CYTHON = False
    fname = "extinction.c"

# Skip importing numpy & cython if we're just doing setup.py egg_info.
if (any('--' + opt in sys.argv for opt in Distribution.display_option_names +
        ['help-commands', 'help']) or len(sys.argv) == 1
    or sys.argv[1] in ('egg_info', 'clean', 'help')):
    extensions=[]
else:
    try:
        import numpy
    except ImportError:
        raise SystemExit("NumPy is required for '{}'".format(sys.argv[1]))

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

    if USE_CYTHON:
        from Cython.Build import cythonize
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
