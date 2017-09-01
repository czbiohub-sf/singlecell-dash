#!/usr/bin/env python

from setuptools import setup

version = '1.0.0'

required = open('requirements.txt').read().split('\n')

setup(
    name='singlecell-dash',
    version=version,
    description='Dashboard for visualizing per-plate sequencing quality metrics and basic',
    author='czbiohub',
    author_email='olga.botvinnik@czbiohub.org',
    url='https://github.com/czbiohub/singlecell-dash',
    packages=['singlecell_dash'],
    install_requires=required,
    long_description='See ' + 'https://github.com/czbiohub/singlecell-dash',
    license='MIT'
)
