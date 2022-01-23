# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>
.. codeauthor:: Zander Harteveld <zandermilanh@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import os
import shutil

# External Libraries
from libconfig import Config
core = Config()

# This Library

with core.ifndef():
    # Register IO control options
    core.register_option('system', 'verbose', False, 'bool', 'Makes topobuilder chatty.')
    core.register_option('system', 'debug', False, 'bool', 'Makes topobuilder VERY chatty.')
    core.register_option('system', 'strict', False, 'bool', 'When True, warings become exit points.')
    core.register_option('system', 'overwrite', False, 'bool', 'Overwrite existing structure files.')
    core.register_option('system', 'forced', False, 'bool', 'Ignore checkpoints and redo calculations.')
    core.register_option('system', 'image', '.png', 'string', 'Format to output images', ['.png', '.svg'])
    core.register_option('system', 'pymol', False, 'bool', 'When True, create pymol scripts to visualize steps.')
    core.register_option('system', 'jupyter', 'JPY_PARENT_PID' in os.environ, 'bool', 'Is TopoBuilder run from a notebook?', locked=True)

    # For plugins that require SLURM submission
    core.register_option('slurm', 'use', True, 'bool', 'Use SLURM cluster submission system.')
    core.register_option('slurm', 'nodes', 1, 'int', 'Nodes to use per job.')
    core.register_option('slurm', 'ntasks-per-node', 1, 'int', 'Tasks to perform in a node.')
    core.register_option('slurm', 'cpus-per-task', 1, 'int', 'CPUs to use per task.')
    core.register_option('slurm', 'memory', 4096, 'int', 'Memory required per node.')
    core.register_option('slurm', 'partition', 'serial', 'string', 'Name of the available SLURM partition.')
    core.register_option('slurm', 'array', 700, 'int', 'Into how may nodes is the search splitted.')
    core.register_option('slurm', 'time', '10:00:00', 'string', 'Expected running time.')
    core.register_option('slurm', 'logs', os.getcwd(), 'path_in', 'Path on were to dump the log files.')

    # For plugins that require MASTER
    core.register_option('master', 'master', shutil.which('master'), 'path_in', 'MASTER executable.')
    core.register_option('master', 'create', shutil.which('createPDS'), 'path_in', 'createPDS executable.')
    core.register_option('master', 'pds', None, 'path_in', 'Local PDS database.')
    core.register_option('master', 'pdb', None, 'path_in', 'Local PDB database.')

    # For plugins that requires RosettaScripts
    core.register_option('rosetta', 'scripts', None, 'path_in', 'Full path to the rosetta_scripts executable.')
    core.register_option('psipred', 'script', None, 'path_in',
                         'Full path to the PSIPRED executable. If present adds PSIPRED to Rosetta reports.')

    # There are different levels of configuration files that can be picked.
    # If any configuration file is set up, the priority goes as follows:
    #   1) Local config file (in the actual executable directory)
    #   2) Root of the current working repository (if any)
    #   3) User's home path
    config_file = core.get_local_config_file('.topobuilder.cfg')
    if config_file is not None:
        core.set_options_from_YAML( config_file )

    core.lock_configuration()
