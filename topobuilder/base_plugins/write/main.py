# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import List, Optional
import sys
from pathlib import Path
from collections import OrderedDict

# External Libraries
import yaml
try:
    from yaml import CDumper as Dumper
except ImportError:
    from yaml import Dumper
from yaml.representer import SafeRepresenter

# This Library
from topobuilder.case import Case
import topobuilder.core as TBcore

__all__ = ['apply']

_PLT_TYPES_ = ['sketchXZ', 'sketchXY']


def apply( cases: List[Case],
           prtid: int,
           prefix: Optional[str] = 'checkpoint',
           **kwargs ) -> List[Case]:
    """Write the Case with all the applied modifications.
    """
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('--- TB PLUGIN: WRITE ---\n')

    counter = '' if prtid is -1 else '.{:02d}'.format(prtid)
    paths = [c.main_path for c in cases]
    try:
        path = paths[0].resolve()
        for p in paths:
            path = path.relative_to(p.resolve())
    except ValueError:
        path = Path.cwd()

    YAML_Dumper()
    ofile = str(path.joinpath(prefix)) + counter + '.yml'
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('Writing checkpoint file {} for {} case(s).\n'.format(ofile, len(cases)))
    with open(ofile, "w") as stream:
        yaml.dump_all([c.data for c in cases], stream, Dumper=Dumper, default_flow_style=False)

    return cases


def YAML_Dumper():
    # This is required for YAML to properly print the Schema as an OrderedDict
    # Adapted from https://gist.github.com/oglops/c70fb69eef42d40bed06 to py3
    def dict_representer(dumper, data):
        return dumper.represent_dict(data.items())

    Dumper.add_representer(OrderedDict, dict_representer)
    Dumper.add_representer(str, SafeRepresenter.represent_str)
