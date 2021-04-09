#!/usr/bin/env python3
import os
from distutils.core import setup
from setuptools import find_packages

author = 'James Robertson'

classifiers = """
Development Status :: 3 - Alpha
Environment :: Console
License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)
Intended Audience :: Science/Research
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Bio-Informatics
Programming Language :: Python
Programming Language :: Python :: 3.7
Programming Language :: Python :: 3.8
Programming Language :: Python :: Implementation :: CPython
Operating System :: POSIX :: Linux
""".strip().split('\n')


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


exec(open('sars_cov_2_kmer_typing/version.py').read())

setup(
    name='sars_cov_2_kmer_typing',
    include_package_data=True,
    version='0.0.1',
    python_requires='>=3.7.0,<4',
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    packages=find_packages(exclude=['tests', 'databases']),
    url='https://github.com/phac-nml/mob-suite',
    license='GPLv3',
    author='James Robertson',
    author_email='james.robertson@canada.ca',
    description=(
        'MOB-suite is a set of tools for finding, typing and reconstruction of plasmids from draft and complete genome assemblies.'),
    keywords='Plasmids finding typing reconstruction',
    classifiers=classifiers,
    package_dir={'sars_cov_2_kmer_typing': 'sars_cov_2_kmer_typing'},
    package_data={'sars_cov_2_kmer_typing': ['config.json']},

    install_requires=[
        'numpy',
        'pandas',
        'biopython',
        'scipy',
        'six',
    ],

    entry_points={
        'console_scripts': [
            'mob_init=mob_suite.mob_init:main',
            'mob_recon=mob_suite.mob_recon:main',
            'mob_cluster=mob_suite.mob_cluster:main',
            'mob_typer=mob_suite.mob_typer:main',
        ],
    },
)
