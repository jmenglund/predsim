#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup
from os.path import join, dirname
from io import open


setup(
    name='predsim',
    version='0.1.0',
    description=(
        'Command-line tool for simulating predictive '
        'datasets from MrBayes\' output.'),
    long_description=open(
        join(dirname(__file__), 'README.rst'), encoding='utf-8').read(),
    py_modules=['predsim'],
    install_requires=['dendropy>=4.0', 'pandas>=0.16'],
    entry_points={
        'console_scripts': [
            'predsim = predsim:main']},
    author='Markus Englund',
    author_email='jan.markus.englund@gmail.com',
    url='https://github.com/jmenglund/predsim',
    license='MIT',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    keywords=['simulation', 'Seq-Gen', 'DendroPy'])
