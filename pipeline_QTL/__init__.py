# Python 3.3 changed the use of __init__.py
# To avoid namespace problems use this for 2 and 3 compatibility with pkgutil
# Otherwise omit __init__.py entirely unless you have sub-packages

# See:
# https://pymotw.com/2/pkgutil/
# https://packaging.python.org/guides/packaging-namespace-packages/

# and examples in:
# https://github.com/pypa/sample-namespace-packages/tree/master/pkgutil

__path__ = __import__('pkgutil').extend_path(__path__, __name__)

# If you have modules in the same directory (as opposed to sub-directories),
# add:
name = 'pipeline_QTL'

# For each sub-directory, add an __init__.py that only contains:
# name = 'sub-directory'
