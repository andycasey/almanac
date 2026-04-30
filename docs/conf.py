# Configuration file for the Sphinx documentation builder.
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here.
sys.path.insert(0, os.path.abspath('../src'))

# -- Project information -----------------------------------------------------

project = 'almanac'
copyright = '2024, Andy Casey'
author = 'Andy Casey'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
try:
    import almanac
    version = almanac.__version__
    release = version
except ImportError:
    version = '0.1.11'
    release = version

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.doctest',
    'sphinx.ext.intersphinx',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
    'sphinx_autodoc_typehints',
    'sphinx_copybutton',
    'myst_parser',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'README.md']

# The suffix(es) of source filenames.
source_suffix = ['.rst', '.md']

# The master toctree document.
master_doc = 'index'

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.
html_theme = 'furo'

# Theme options are theme-specific and customize the look and feel of a theme
# further.
html_theme_options = {
    "source_repository": "https://github.com/andycasey/almanac/",
    "source_branch": "main",
    "source_directory": "docs/",
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']


# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'almanacdoc'

# -- Extension configuration -------------------------------------------------

# -- Options for intersphinx extension ---------------------------------------

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'astropy': ('https://docs.astropy.org/en/stable/', None),
    'click': ('https://click.palletsprojects.com/en/8.1.x/', None),
}

# -- Options for autodoc extension -------------------------------------------

# This value controls how to represent typehints.
autodoc_typehints = 'description'
autodoc_typehints_description_target = 'documented'

# Mock imports for packages that are not available during build
autodoc_mock_imports = [
    'sdssdb',
    'pydl',
    'colorlog',
    'rich',
    'pydantic',
    'pydantic_core'
]

# -- Options for Napoleon extension ------------------------------------------

napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True

# -- Options for MyST parser ------------------------------------------------

myst_enable_extensions = [
    'deflist',
    'tasklist',
    'html_admonition',
    'html_image',
    'colon_fence',
    'smartquotes',
    'replacements',
    'linkify',
    'strikethrough',
    'substitution',
]

# -- Options for copybutton extension ---------------------------------------

copybutton_prompt_text = r'>>> |\.\.\. |\$ |In \[\d*\]: | {2,5}\.\.\.: | {5,8}: '
copybutton_prompt_is_regexp = True