#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import os
import numpy
import nibabel

# Global parameters
POSSIBLE_AXES_ORIENTATIONS = [
    "LAI", "LIA", "ALI", "AIL", "ILA", "IAL",
    "LAS", "LSA", "ALS", "ASL", "SLA", "SAL",
    "LPI", "LIP", "PLI", "PIL", "ILP", "IPL",
    "LPS", "LSP", "PLS", "PSL", "SLP", "SPL",
    "RAI", "RIA", "ARI", "AIR", "IRA", "IAR",
    "RAS", "RSA", "ARS", "ASR", "SRA", "SAR",
    "RPI", "RIP", "PRI", "PIR", "IRP", "IPR",
    "RPS", "RSP", "PRS", "PSR", "SRP", "SPR"
]
CORRECTION_MATRIX_COLUMNS = {
    "R": (1, 0, 0),
    "L": (-1, 0, 0),
    "A": (0, 1, 0),
    "P": (0, -1, 0),
    "S": (0, 0, 1),
    "I": (0, 0, -1)
}


def swap_affine(axes):
    """ Build a correction matrix, from the given orientation of axes to RAS.

    Parameters
    ----------
    axes: str (manadtory)
        the given orientation of the axes.

    Returns
    -------
    rotation: array (4, 4)
        the correction matrix.
    """
    rotation = numpy.eye(4)
    rotation[:3, 0] = CORRECTION_MATRIX_COLUMNS[axes[0]]
    rotation[:3, 1] = CORRECTION_MATRIX_COLUMNS[axes[1]]
    rotation[:3, 2] = CORRECTION_MATRIX_COLUMNS[axes[2]]
    return rotation


def reorient_image(in_file, axes="RAS", prefix="swap", outdir=None):
    """ Rectify the orientation of an image.

    Parameters
    ----------
    in_file: str (mandatory)
        the input image.
    axes: str (optional, default 'RAS')
        orientation of the original axes X, Y, and Z
        specified with the following convention: L=Left, R=Right,
        A=Anterion, P=Posterior, I=Inferior, S=Superior.
    prefix: str (optional, default 'swap')
        prefix of the output image.
    outdir: str (optional, default None)
        the output directory where the rectified image is saved.
        If None use the same directory as the input image.

    Returns
    -------
    out_file: str
        the rectified image.

    Examples
    --------

    >>> from monkeyfmri.utils.reorientation import reorient_image
    >>> rectified_image = reorient_image('image.nii', 'RAS', 's', None)
    """
    # Check the input image exists on the file system
    if not os.path.isfile(in_file):
        raise ValueError("'{0}' is not a valid filename.".format(in_file))

    # Check that the outdir is valid
    if outdir is not None:
        if not os.path.isdir(outdir):
            raise ValueError("'{0}' is not a valid directory.".format(outdir))
    else:
        outdir = os.path.dirname(in_file)

    # Check that a valid input axes is specified
    if axes not in POSSIBLE_AXES_ORIENTATIONS:
        raise ValueError("Wrong coordinate system: {0}.".format(axes))

    # Get the transformation to the RAS space
    rotation = swap_affine(axes)
    det = numpy.linalg.det(rotation)
    if det != 1:
        raise Exception("Rotation matrix determinant must be one "
                        "not '{0}'.".format(det))

    # Load the image to rectify
    image = nibabel.load(in_file)

    # Get the input image affine transform
    affine = image.get_affine()

    # Apply the rotation to set the image in the RAS coordiante system
    transformation = numpy.dot(rotation, affine)
    image.set_qform(transformation)
    image.set_sform(transformation)

    # Save the rectified image
    fsplit = os.path.split(in_file)
    out_file = os.path.join(outdir, prefix + fsplit[1])
    nibabel.save(image, out_file)

    return out_file
