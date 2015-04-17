#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# CAPSUL import
try:
    from capsul.utils.pilot import pilotfunction
except:
    def pilotfunction(func):
        return func


@pilotfunction
def pilot_qa():
    """ 
    Quanlity assurance
    ===================
    """
    # Pilot imports
    import os
    from caps.toy_datasets import get_sample_data
    from capsul.study_config import StudyConfig
    from capsul.process import get_process_instance

    """
    Study configuration
    -------------------

    We first define the working directory and guarantee this folder exists on
    the file system:
    """
    working_dir = "/volatile/nsap/pclinfmri/qa"
    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)

    """
    And then define the study configuration (here we activate the smart
    caching module that will be able to remember which process has already been
    processed):
    """
    study_config = StudyConfig(
        modules=["SmartCachingConfig"],
        use_smart_caching=True,
        output_directory=working_dir)

    """
    Load the toy dataset
    --------------------

    To do so, we use the get_sample_data function to download the toy
    dataset on the local file system (here localizer data):
    """
    toy_dataset = get_sample_data("localizer")

    """
    The toy_dataset is an Enum structure with some specific elements of
    interest:

        * **??**: ??.

    Processing definition
    ---------------------

    First create the
    :ref:`slice timing pipeline <pclinfmri.preproc.pipeline.SliceTiming>` that
    define the different step of the processings:
    """
    pipeline = get_process_instance(
        "pclinfmri.quality_assurance.fmri_quality_assurance.xml")
    print pipeline.get_input_spec()

    """
    Now we need now to parametrize this pipeline:
    """
    pipeline.image_file = toy_dataset.fmri
    pipeline.repetition_time = toy_dataset.TR

    """
    The pipeline is now ready to be run:
    """
    study_config.run(pipeline, executer_qc_nodes=True, verbose=1)

    """
    Results
    -------

    Finally, we print the pipeline outputs:
    """
    print("\nOUTPUTS\n")
    for trait_name, trait_value in pipeline.get_outputs().items():
        print("{0}: {1}".format(trait_name, trait_value))


if __name__ == "__main__":
    pilot_qa()


