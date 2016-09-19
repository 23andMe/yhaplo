#!/usr/bin/env python
import sys

from setuptools import setup, find_packages

PY2 = sys.version_info[0] == 2
PY3 = sys.version_info[0] == 3

if PY2:
    readme = open('README.md').read()
    license = ''  # not yet: open('LICENSE').read()
elif PY3:
    readme = open('README.md', encoding='utf-8').read()
    license = ''  # not yet: open('LICENSE', encoding='utf-8').read()

setup(
    name='yhaplo',
    version='1.0.0',
    url='https://github.com/23andMe/yhaplo',
    author='David Poznik',
    author_email='dpoznik@23andme.com',
    description='Y Chromosome Haplogroup Calling',
    long_description=readme,
    license='TBD',  # TODO
    packages=find_packages(),
    include_package_data=True,
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'License :: Free for non-commercial use',  # TODO ?
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
    ]
)
