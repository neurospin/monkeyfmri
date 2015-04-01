#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# CASPUL import
from capsul.utils.function_to_process import register_processes

# MONKEYFMRI import
from .reorientation import reorient_image


# Register new processes
register_processes([reorient_image])
