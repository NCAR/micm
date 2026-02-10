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
import subprocess
sys.path.insert(0, os.path.abspath('.'))

DOCS_SOURCE_DIR = os.path.abspath(os.path.dirname(__file__))
REPO_ROOT_DIR = os.path.abspath(os.path.join(DOCS_SOURCE_DIR, '..', '..'))
BUILD_DIR = os.path.join(REPO_ROOT_DIR, 'build')
DOXYGEN_XML_DIR = os.path.join(BUILD_DIR, 'docs', 'doxygen', 'xml')


def _run_command(command, cwd=None):
    try:
        subprocess.run(command, cwd=cwd, check=True)
    except (OSError, subprocess.CalledProcessError) as exc:
        sys.stderr.write(f"Command failed: {command}\n{exc}\n")
        raise


def _ensure_doxygen_xml():
    read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'
    if not read_the_docs_build:
        return

    cache_file = os.path.join(BUILD_DIR, 'CMakeCache.txt')
    if not os.path.exists(cache_file):
        _run_command([
            'cmake',
            '-S', REPO_ROOT_DIR,
            '-B', BUILD_DIR,
            '-D', 'MICM_BUILD_DOCS=ON'
        ])

    _run_command([
        'cmake',
        '--build', BUILD_DIR,
        '--target', 'Doxygen'
    ])

# -- Project information -----------------------------------------------------


project = 'MICM'
copyright = f'2022-{datetime.datetime.now().year}, NSF-NCAR/ACOM'
author = 'NSF-NCAR/ACOM'

regex = r'project\(\w+\s+VERSION\s+(\d+\.\d+\.\d+)'
version = '0.0.0'
# read the version from the cmake files
with open(f'../../CMakeLists.txt', 'r') as f:
    for line in f:
        match = re.match(regex, line)
        if match:
            version = match.group(1)
release = f'{version}'

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
breathe_projects = {
    "micm": DOXYGEN_XML_DIR
}

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
    'mb': ('https://ncar.github.io/music-box/', None)
}

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'pydata_sphinx_theme'

html_theme_options = {
    "external_links": [],
    "github_url": "https://github.com/NCAR/micm",
    "navbar_end": ["navbar-icon-links"],
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


def setup(app):
    app.connect("builder-inited", lambda _app: _ensure_doxygen_xml())
