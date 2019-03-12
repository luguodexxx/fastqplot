#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/3/5 5:20 PM
__author__ = 'Zhou Ran'

from setuptools import setup, find_packages
from fastqs.version import __version__

setup(name='fastqs',
      version=__version__,
      python_requires='>3.6',
      description='fastqs',
      author='Ran zhou',
      author_email='ranzhou1005@gmail.com',
      url='https://github.com/luguodexxx/fastqs',
      license='MIT',
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3',
      ],
      keywords='sashimiplot',
      packages=find_packages(),
      install_requires=[
          'requests',
          'click',
          'matplotlib'
      ],
      entry_points={
          'console_scripts': [
              'fqplot=fastqs.cli:cli'
          ],
      },
      )
