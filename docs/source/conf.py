# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import datetime
import re
sys.path.insert(0, os.path.abspath('.'))

# -- Project information -----------------------------------------------------

project = 'MICM'
copyright = f'2022-{datetime.datetime.now().year}, NCAR/UCAR'
author = 'NCAR/UCAR'

suffix = os.getenv("SWITCHER_SUFFIX", "")
# the suffix is required. This is controlled by the dockerfile that builds the docs
regex = r'project\(\w+\s+VERSION\s+(\d+\.\d+\.\d+)'
version = '0.0.0'
# read the version from the cmake files
with open(f'../../CMakeLists.txt', 'r') as f:
    for line in f:
        match = re.match(regex, line)
        if match:
            version = match.group(1)
release = f'{version}{suffix}'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
  'breathe',
  'sphinx_copybutton',
  'sphinx_design',
  'sphinxcontrib.bibtex',
  'sphinx.ext.intersphinx'
]

bibtex_bibfiles = ['references.bib']
suppress_warnings = ["bibtex.missing_field"]

breathe_default_project = "micm"

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# -- Link to Intersphinx Mappings ---

intersphinx_mapping = {
    'musica': ('https://ncar.github.io/musica/', None),
    'mc': ('https://ncar.github.io/MechanismConfiguration/', None),
    'mb': ('https://ncar.github.io/music-box/',None)
}

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'pydata_sphinx_theme'

html_theme_options = {
    "external_links": [],
    "github_url": "https://github.com/NCAR/micm",
    "navbar_end": ["version-switcher", "navbar-icon-links"],
    "switcher": {
        "json_url": "https://ncar.github.io/micm/_static/switcher.json",
        "version_match": release,
    },
   "pygment_light_style": "tango",
   "pygment_dark_style": "monokai"
}

html_css_files = [
    'custom.css'
]

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_favicon = '_static/favicon/favicon.ico'