
import os
import sys
sys.path.insert(0, os.path.abspath('.'))

project = 'MultiDIC'
author = 'UCLouvain - INMA - COSY'
release = '1.0'

extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon', 'sphinx.ext.viewcode']
templates_path = ['_templates']
exclude_patterns = []
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
