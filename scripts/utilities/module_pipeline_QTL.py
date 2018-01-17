'''
Utilities for |project_name|
============================

:Authors: |author_names|
:Release: |version|
:Date: |today|


Purpose
=======

|description|


Documentation
=============

    For more information see:

        |url|

'''
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import

from future import standard_library
standard_library.install_aliases()
import os
import sys

# Set up calling parameters from INI file:
# Modules with Py2 to 3 conflicts
try:
    import configparser
except ImportError:  # Py2 to Py3
    import ConfigParser as configparser
# Global variable for configuration file ('.ini'):
CONFIG = configparser.ConfigParser(allow_no_value = True)
#################


################# 
cwd = os.getcwd()
def getINIdir(path = cwd):
    ''' Search for an INI file, default is where the current working directory '''
    f_count = 0
    for f in os.listdir(path):
        if (f.endswith('.ini') and not f.startswith('tox')):
            f_count += 1
            INI_file = f
    if f_count == 1:
        INI_file = os.path.abspath(os.path.join(path, INI_file))
    elif (f_count > 1 or f_count == 0):
        INI_file = os.path.abspath(path)
        print('You have no project configuration (".ini") file or more than one',
              'in the directory:', '\n', path)

    return(INI_file)
################# 


#################
# Get locations of source code
    # os.path.join note: a subsequent argument with an '/' discards anything
    # before it
    # For function to search path see: 
    # http://stackoverflow.com/questions/4519127/setuptools-package-data-folder-location
# MANIFEST.in file instructs the project_quickstart/templates folder to be included in installation

_ROOT = os.path.abspath(os.path.dirname(__file__))
def getDir(path = _ROOT):
    ''' Get the absolute path to where this function resides. Useful for
    determining the user's path to a package. If a sub-directory is given it
    will be added to the path returned. Use '..' to go up directory levels. '''
   # src_top_dir = os.path.abspath(os.path.join(_ROOT, '..'))
    src_dir = _ROOT
    return(os.path.abspath(os.path.join(src_dir, path)))
#################
