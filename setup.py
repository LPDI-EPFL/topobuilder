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
    url='https://github.com/jaumebonet/RosettaSilentToolbox',

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
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Environment :: Console',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Operating System :: MacOS',
        'Operating System :: Unix',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],

    project_urls={
        'Source': 'https://github.com/lpdi-epfl/topobuilder',
        'Tracker': 'https://github.com/lpdi-epfl/topobuilder/issues',
    },

    platforms='UNIX',
    keywords='development',

    install_requires=[x.strip() for x in open('REQUIREMENTS').readlines()],

    packages=find_packages(exclude=['docs', 'demo', 'sphinx-docs']),
    include_package_data=True,
    package_data={
        'topobuilder': ['REQUIREMENTS', ]
    },
    cmdclass=versioneer.get_cmdclass()
)
