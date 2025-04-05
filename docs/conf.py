
import os
import sys
sys.path.insert(0, os.path.abspath('.'))

project = 'MultiDIC'
author = 'UCLouvain - INMA - COSY'
release = '1.0'

extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon', 'sphinx.ext.viewcode']
templates_path = ['_templates']
exclude_patterns = []
import sphinx_rtd_theme

html_theme = 'sphinx_rtd_theme'
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
html_static_path = ['_static']
