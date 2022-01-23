# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>
.. codeauthor:: Zander Harteveld <zandermilanh@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries

# External Libraries

# This Library
from topobuilder.core import core

core.open = True  # Open config to add more options.
with core.ifndef():
    # For plugins Analytics
    core.register_option('statistics', 'molprobity',  None, 'path_in',
                         'Full path to the MolProbity oneline-analysis executable. If present adds it to statistics.')
    core.register_option('statistics', 'tmalign',  None, 'path_in', 'Full path to TMalign executable.')
    core.register_option('statistics', 'trrosetta_repo', None, 'path_in',
                         'Full path to the trRosetta folder. If present adds it to statistics.')
    core.register_option('statistics', 'trrosetta_wts', None, 'path_in', 'Full path to the trRosetta pretrained weights.')
    core.register_option('statistics', 'trrosetta_env', None, 'path_in', 'Full path to conda environment to use for trRosetta.')

    # There are different levels of configuration files that can be picked.
    # If any configuration file is set up, the priority goes as follows:
    #   1) Local config file (in the actual executable directory)
    #   2) Root of the current working repository (if any)
    #   3) User's home path
    config_file = core.get_local_config_file('.topobuilder.cfg')
    if config_file is not None:
        core.set_options_from_YAML( config_file )

    core.lock_configuration()
