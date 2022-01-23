import os
from setuptools import setup, find_packages
import versioneer


def read_file(path):
    with open(os.path.join(os.path.dirname(__file__), path)) as fp:
        return fp.read()


setup(
    name='topobuilder',
    version=versioneer.get_version(),

    description='Automatic Building of Functional de novo Topologies',
    long_description=read_file('README.rst'),

    # The project's main homepage.
    url='https://github.com/jaumebonet/topobuilder',

    # Author details
    author='Jaume Bonet',
    author_email='jaume.bonet@gmail.com',

    # Choose your license
    license='MIT',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Environment :: Console',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Operating System :: MacOS',
        'Operating System :: Unix',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],

    project_urls={
        'Documentation': 'http://jaumebonet.cat/topobuilder',
        'Source': 'https://github.com/jaumebonet/topobuilder/',
        'Tracker': 'https://github.com/jaumebonet/topobuilder/issues',
    },

    platforms='UNIX',
    keywords='development',

    install_requires=[x.strip() for x in open('REQUIREMENTS').readlines()],

    packages=find_packages(exclude=['docs', 'demo']),
    include_package_data=True,
    package_data={
        'topobuilder': ['REQUIREMENTS', ]
    },
    entry_points={
        'console_scripts':
            ['topo.case=topobuilder.interface.cli.case:cli_case_template',
             'topo.absolute=topobuilder.interface.cli.case:cli_absolute_case',
             'topo.builder=topobuilder.interface.cli.actions:cli_build',
             'topo.protocol=topobuilder.interface.cli.actions:cli_protocol',
             ]
    },
    cmdclass=versioneer.get_cmdclass(),
)
