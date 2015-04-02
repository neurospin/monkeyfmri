#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################


"""SPM wrappers for preprocessing data
"""

__docformat__ = 'restructuredtext'

# Standard library imports
from copy import deepcopy
import os

# Third-party imports
import numpy as np

# Local imports
from nipype.interfaces.base import (OutputMultiPath, TraitedSpec, isdefined,
                                    traits, InputMultiPath, File)
from nipype.interfaces.spm.base import (SPMCommand, scans_for_fname,
                                        func_is_3d,
                                        scans_for_fnames, SPMCommandInputSpec)
from nipype.utils.filemanip import (fname_presuffix, filename_to_list,
                                    list_to_filename, split_filename)


class ApplyDeformationFieldInputSpec(SPMCommandInputSpec):
    """
    Parameters
    ----------
    in_files : list of str (mandatory)
        Files on which the deformation is applied.

    deformation_field : str (mandatory)
        SN SPM deformation file.

    bounding_box : list of list of float
        3x2-element list of lists (opt).

    voxel_sizes : list of float
        3-element list (opt).

    interpolation : int
        Degree of b-spline used for interpolation (from 0 to 7).
    """
    in_files = InputMultiPath(
        File(exists=True),
        mandatory=True,
        field='fnames',
        desc='Files on which deformation is applied')
    deformation_field = File(
        exists=True,
        mandatory=True,
        field='comp{1}.def',
        desc='SN SPM deformation file')
    bounding_box = traits.List(
        traits.List(traits.Float(),
        minlen=3, maxlen=3),
        field='comp{2}.idbbvox.bb',
        minlen=2, maxlen=2,
        desc='3x2-element list of lists (opt)')
    voxel_sizes = traits.List(
        traits.Float(),
        [1., 1., 1.],
        field='comp{2}.idbbvox.vox',
        minlen=3, maxlen=3,
        desc='3-element list (opt)')
    interpolation = traits.Range(
        low=0,
        high=7,
        field='interp',
        desc='degree of b-spline used for interpolation')
    out_prefix = traits.String(
        'w', field='prefix',
        usedefault=True,
        desc='aplly deformation field output prefix')


class ApplyDeformationFieldOutputSpec(TraitedSpec):
    """
    Returns
    -------
    normalized_files : list of str
        Converted files.
    """
    normalized_files = OutputMultiPath(
        File(exists=True),
        desc='converted files')


class ApplyDeformationField(SPMCommand):
    """ Uses SPM to apply inverse deformation field to given files.

    Examples
    --------
    >>> import nsap.nipype.spm_interfaces as spm
    >>> f = spm.ApplyformationField()
    >>> f.inputs.in_files = 'functional.nii'
    >>> f.inputs.deformation_field = 'y_t1_localizer.nii'
    >>> f.run()
    """

    input_spec = ApplyDeformationFieldInputSpec
    output_spec = ApplyDeformationFieldOutputSpec

    _jobtype = 'util'
    _jobname = 'defs'

    def _format_arg(self, opt, spec, val):
        """Convert input to appropriate format for SPM
        """
        if opt in ['in_files']:
            return scans_for_fnames(filename_to_list(val), keep4d=False)
        if opt == 'deformation_field':
            return np.array([list_to_filename(val)], dtype=object)
        return super(ApplyDeformationField, self)._format_arg(opt, spec, val)

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['normalized_files'] = []
        for filename in self.inputs.in_files:
            _, fname = os.path.split(filename)
            outputs['normalized_files'].append(
                os.path.join(self.inputs.output_directory,
                             '%s%s' % (self.inputs.out_prefix, fname)))
        return outputs
