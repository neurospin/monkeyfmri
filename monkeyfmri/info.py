##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

""" This file contains parameters for the project.
"""

_version_major = 0
_version_minor = 0
_version_micro = 1

# Format expected by setup.py and doc/source/conf.py: string of form "X.Y.Z"
__version__ = "{0}.{1}.{2}".format(_version_major, _version_minor,
                                   _version_micro)

CLASSIFIERS = ["Development Status :: 3 - Alpha",
               "Environment :: Console",
               "Operating System :: OS Independent",
               "Programming Language :: Python",
               "Topic :: Scientific/Engineering",
               "Topic :: Utilities"]

description = "Monkey FMRI"

long_description = """\
Monkey FMRI
===========


ToDo
"""

# versions for dependencies
NIBABEL_MIN_VERSION = "1.3.0"
NUMPY_MIN_VERSION = "1.3"

# Main setup parameters
NAME = "monkeyfmri"
MAINTAINER = "Bechir Jarraya"
MAINTAINER_EMAIL = "bechir.jarraya@gmail.com"
DESCRIPTION = description
LONG_DESCRIPTION = long_description
URL = "https://github.com/neurospin/monkeyfmri"
DOWNLOAD_URL = ""
LICENSE = "CeCILL-B"
CLASSIFIERS = CLASSIFIERS
AUTHOR = "monkeyfmri developers"
AUTHOR_EMAIL = "bechir.jarraya@gmail.com"
PLATFORMS = "OS Independent"
MAJOR = _version_major
MINOR = _version_minor
MICRO = _version_micro
VERSION = __version__
PROVIDES = ["monkeyfmri"]
REQUIRES = ["numpy>={0}".format(NUMPY_MIN_VERSION),
            "nibabel>={0}".format(NIBABEL_MIN_VERSION)]
EXTRA_REQUIRES = {"doc": ["sphinx>=1.0"]}
