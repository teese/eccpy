from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'readme.rst')) as f:
    long_description = f.read()

classifiers = """\
Intended Audience :: Science/Research
License :: OSI Approved :: MIT License
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
      license='MIT',
      packages=find_packages(),
      classifiers=classifiers.splitlines(),
      install_requires=["pandas", "numpy", "matplotlib", "pytest", "scipy", "xlrd", "openpyxl", "pytest"],
      keywords="EC50 LD50 IC50 doseresponse concentration dose inhibitor sigmoidal curve",
      # obtains package data from MANIFEST.in
      include_package_data=True,
      version="0.5.0")
