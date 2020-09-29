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
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'readme.rst')) as f:
    long_description = f.read()

classifiers = """\
Intended Audience :: Science/Research
License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)
Programming Language :: Python :: 3
Topic :: Scientific/Engineering :: Bio-Informatics
Topic :: Scientific/Engineering :: Medical Science Apps.
Topic :: Scientific/Engineering :: Chemistry
"""

setup(name='eccpy',
      description="High-throughput calculation of EC50 values.",
      author="Mark Teese",
      author_email="mark.teese@checkthereadme.com",
      url="https://github.com/teese/eccpy",
      download_url='https://github.com/teese/eccpy/archive/0.5.0.tar.gz',
      project_urls={'Wiki': 'https://github.com/teese/eccpy/wiki',
                    'LangoschLab': 'http://cbp.wzw.tum.de/index.php?id=9',
                    "TU_Muenchen": "https://www.tum.de",
                    "TNG Technology Consulting GmbH": "https://www.tngtech.com/en/index.html"},
      long_description=long_description,
      long_description_content_type='text/x-rst',
      license='LGPLv3',
      packages=find_packages(),
      classifiers=classifiers.splitlines(),
      install_requires=["pandas", "numpy", "matplotlib", "pytest", "scipy", "xlrd", "openpyxl", "pytest"],
      keywords="EC50 LD50 IC50 doseresponse concentration dose inhibitor sigmoidal curve",
      # obtains package data from MANIFEST.in
      include_package_data=True,
      version="0.5.0")
