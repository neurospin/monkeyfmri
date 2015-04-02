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
from .io import element_to_list
from .io import list_to_element
from .io import ungzip_file
from .io import gzip_file
from .io import spm_tissue_probability_maps

# Register new processes
register_processes([element_to_list, list_to_element, ungzip_file, gzip_file,
                    spm_tissue_probability_maps])
