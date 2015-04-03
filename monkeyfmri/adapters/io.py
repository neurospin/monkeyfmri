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
import gzip


def element_to_list(element):
    """ Set an element to an empty list.

    <process>
        <return name="adaptedelement" type="List_Any" desc="the returned
            list containing the input element only."/>
        <input name="element" type="Any" desc="an input element."/>
    </process>
    """
    return [element]


def list_to_element(listobj):
    """ Get the singleton list element.

    <process>
        <return name="element" type="Any" desc="the returned
            list single element."/>
        <input name="listobj" type="List_Any" desc="an input singleton list."/>
    </process>
    """
    # Check that we have a singleton list
    if len(listobj) != 1:
        raise ValueError("A list with '{0}' element(s) is not a singleton "
                         "list.".format(len(listobj)))
    return listobj[0]


def ungzip_file(fname, prefix="u", output_directory=None):
    """ Copy and ungzip the input file.

    <process>
        <return name="ungzipfname" type="File" desc="the returned
            ungzip file."/>
        <input name="fname" type="File" desc="an input file to ungzip."/>
        <input name="prefix" type="String" desc="the prefix of the result
            file."/>
        <input name="output_directory" type="Directory" desc="the output
            directory where ungzip file is saved."/>
    </process>
    """
    # Check the input file exists on the file system
    if not os.path.isfile(fname):
        raise ValueError("'{0}' is not a valid filename.".format(fname))

    # Check that the outdir is valid
    if output_directory is not None:
        if not os.path.isdir(output_directory):
            raise ValueError(
                "'{0}' is not a valid directory.".format(output_directory))
    else:
        output_directory = os.path.dirname(fname)

    # Get the file descriptors
    base, extension = os.path.splitext(fname)
    basename = os.path.basename(base)

    # Ungzip only known extension
    if extension in [".gz"]:

        # Generate the output file name
        basename = prefix + basename
        ungzipfname = os.path.join(output_directory, basename)

        # Read the input gzip file
        with gzip.open(fname, "rb") as gzfobj:
            data = gzfobj.read()

        # Write the output ungzip file
        with open(ungzipfname, "w") as openfile:
            openfile.write(data)

    # Default, unknown compression extension: the input file is returned
    else:
        ungzipfname = fname

    return ungzipfname


def gzip_file(fname, prefix="g", output_directory=None):
    """ Copy and gzip the input file.

    <process>
        <return name="gzipfname" type="File" desc="the returned
            gzip file."/>
        <input name="fname" type="File" desc="an input file to gzip."/>
        <input name="prefix" type="String" desc="the prefix of the result
            file."/>
        <input name="output_directory" type="Directory" desc="the output
            directory where gzip file is saved."/>
    </process>
    """
    # Check the input file exists on the file system
    if not os.path.isfile(fname):
        raise ValueError("'{0}' is not a valid filename.".format(fname))

    # Check that the outdir is valid
    if output_directory is not None:
        if not os.path.isdir(output_directory):
            raise ValueError(
                "'{0}' is not a valid directory.".format(output_directory))
    else:
        output_directory = os.path.dirname(fname)

    # Get the file descriptors
    base, extension = os.path.splitext(input_image)

    # Ungzip only non compressed file
    if extension not in [".gz"]:

        # Generate the output file name
        basename = prefix + base + ".gz"
        gzipfname = os.path.join(output_directory, basename)

        # Write the output gzip file
        with open(fname, "rb") as openfile:
            with gzip.open(gzipfname, "w") as gzfobj:
                gzfobj.writelines(openfile)

    # Default, the input file is returned
    else:
        gzipfname = fname

    return gzipfname


def spm_tissue_probability_maps():
    """ SPM tissue probability maps.

    <process>
        <return name="tmp_struct" type="List_Any" desc="a struct containing the
            spm tissue probability map descriptions."/>
    </process>
    """
    # Try to import the resource
    try:
        from caps.toy_datasets import get_sample_data
    except:
        raise ImportError("Can't import 'caps'.")
    tmp_file = get_sample_data("tpm").all

    # Format the tpm for spm
    tissue1 = ((tmp_file, 1), 2, (True, True), (False, True))
    tissue2 = ((tmp_file, 2), 2, (True, True), (False, True))
    tissue3 = ((tmp_file, 3), 2, (True, False), (False, False))
    tissue4 = ((tmp_file, 4), 3, (False, False), (False, False))
    tissue5 = ((tmp_file, 5), 4, (False, False), (False, False))

    return [tissue1, tissue2, tissue3, tissue4, tissue5]


