"""setup.py

This file configures the package for installation.
"""

import os
import inspect
import sys
import time

from setuptools import setup, find_packages

# The versioning here follows semantic versioning 2.0.0
# see https://semver.org
MAJOR = 0
MINOR = 0
PATCH = 1

DESCRIPTION = "Functions for generating Dendro C++ Code and Projects"

FILE_DIR = os.path.dirname(__file__)

# get the requirements file
requires_file = os.path.join(FILE_DIR, 'requirements.txt')
with open(requires_file, 'r') as file_handle:
    REQURES = file_handle.readlines()

# pull the requirements out
REQURES = [r.replace('\n', '') for r in REQURES]

with open(os.path.join(FILE_DIR, "README.md"), "r") as file_handle:
    LONG_DESCRIPTION = file_handle.read()

setup(name="dendrosym",
      version=f"{MAJOR}.{MINOR}.{PATCH}",
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      packages=find_packages(exclude=["bssn", "GR", "tests", "data", "docs"]),
      install_requires=REQURES,
      author="",
      python_requires=">=3.6",
      classifiers=["Programming Language :: Python :: 3"])
