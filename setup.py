"""
yhaplo identifies Y-chromosome haplogroups

Sequence data will yield the most highly resolved classifications,
but the algorithm also works well with chip-based genotype data,
provided a reasonable number of phylogenetically informative sites
have been assayed.
"""

import setuptools
import subprocess

DOC_LIST = __doc__.split('\n')

def get_version_txt():
    'reads version from text file'

    from yhaplo import __version__
    return __version__

def get_version_git():
    'extracts version from git tag'

    checksum = subprocess.check_output(
        'git rev-list --tags --max-count=1'.split()).strip().decode('utf-8')
    version = (subprocess.check_output(
        'git describe --tags {}'.format(checksum).split())
        .strip().decode('utf-8'))

    return version

setuptools.setup(
    name='yhaplo',
    version=get_version_txt(),
    author='David Poznik',
    description=DOC_LIST[1],
    long_description='\n'.join(DOC_LIST[3:]),
    license='https://github.com/23andMe/yhaplo/blob/master/LICENSE.txt',
    url='https://github.com/23andMe/yhaplo',
    packages=setuptools.find_packages(),
    include_package_data=True,
    zip_safe=True,
    install_requires=['six>=1.12'],
    setup_requires=[],
    entry_points = {
        'console_scripts': [
            'yhaplo=yhaplo.call_haplogroups:call_haplogroups',
            'convert_to_genos=yhaplo.convert_to_genos:main',
            'plot_tree=yhaplo.plot_tree:main',
        ],
    }
)
