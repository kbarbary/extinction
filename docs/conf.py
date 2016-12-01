# -*- coding: utf-8 -*-
#
# build requires sphinx_rtd_theme and numpydoc.

import sys
import os
import sphinx_rtd_theme
import matplotlib.sphinxext.plot_directive
import extinction

# ensure that plot helper is on the path
sys.path.insert(0, os.path.abspath(__file__))

# generate api directory if it doesn't already exist
if not os.path.exists('api'):
    os.mkdir('api')

# -- General configuration ------------------------------------------------

intersphinx_mapping = {
    'python': ('http://docs.python.org/', None),
    'numpy': ('http://docs.scipy.org/doc/numpy/', None)}

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'numpydoc',
    matplotlib.sphinxext.plot_directive.__name__]

numpydoc_show_class_members = False
autosummary_generate = True
autoclass_content = "class"
autodoc_default_flags = ["members", "no-special-members"]

# The suffix of source filenames.
source_suffix = '.rst'


# The master toctree document.
master_doc = 'index'

# General information about the project.
project = u'extinction'
copyright = u'2016, Kyle Barbary and contributors'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# THe short X.Y version.
version = '.'.join(extinction.__version__.split('.')[0:2])

# The full version, including alpha/beta/rc tags.
release = extinction.__version__

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_build']

# The reST default role (used for this markup: `text`) to use for all
# documents.
default_role = 'obj'

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'


# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = "sphinx_rtd_theme"
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = []

# Output file base name for HTML help builder.
htmlhelp_basename = 'extinctiondoc'


# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
# The paper size ('letterpaper' or 'a4paper').
#'papersize': 'letterpaper',

# The font size ('10pt', '11pt' or '12pt').
#'pointsize': '10pt',

# Additional stuff for the LaTeX preamble.
#'preamble': '',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
  ('index', 'extinction.tex', u'extinction Documentation',
   u'Kyle Barbary', 'manual'),
]

# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    ('index', 'extinction', u'extinction Documentation',
     [u'Kyle Barbary'], 1)
]

# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
  ('index', 'extinction', u'extinction Documentation',
   u'Kyle Barbary', 'extinction', 'One line description of project.',
   'Miscellaneous'),
]
