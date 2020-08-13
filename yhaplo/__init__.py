import os

version_fn = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    'version.txt')

__version__ = open(version_fn).readline().strip()
