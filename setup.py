#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ECCpy is a user-friendly program for the high-throughput calculation of EC50 values.

Copyright (C) 2016  Mark George Teese

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
from setuptools import setup, find_packages

classifiers = """\
Development Status :: Experimental
Intended Audience :: Science/Research
License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)
Programming Language :: Python :: 3
Topic :: Scientific/Engineering :: Bio-Informatics
Topic :: Scientific/Engineering :: Medical Science Apps.
Topic :: Scientific/Engineering :: Chemistry
Classifier: Operating System :: OS Independent
"""

setup(name='eccpy',
      author='Mark Teese',
      author_email='mark.teese /at/ tum.de',
      license='LGPLv3',
      packages=find_packages(),
      classifiers=classifiers.splitlines(),
      platforms=['ALL'],
      keywords=["EC50", "LD50", "IC50", "dose-response", "effective concentration",
                "lethal dose", "inhibitor concentration", "sigmoidal curve"],
      requires=['pandas', 'matplotlib', 'numpy', 'scipy'],
      version='0.3.4')