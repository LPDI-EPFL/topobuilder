# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. func:: build
"""
# Standard Libraries
import argparse
import os
import shutil
from pathlib import Path
from collections import OrderedDict
from typing import Optional, Dict, Union

# External Libraries


# This Library
from topobuilder.case import Case
from topobuilder.io import setup_build
from topobuilder.coordinates import GeneralArchitect

__all__ = ['build']


def build( case: Union[str, Path, Dict],
           overwrite: Optional[bool] = False ):
    """
    """
    info = []
    # 1. Load case and make sure it is absolute.
    if isinstance(case, str):
        case = Path(case)
    data = Case(case).cast_absolute()

    # 2. Create output working directory tree
    paths = setup_build(data, overwrite)

    # 3. Generate Sketch
    arch = GeneralArchitect(data, paths)
    info.append(arch.build_sketch())

    # 4. Calculate Connectivities
    # @TODO

    # 5. Build Connectivities
    info.extend(arch.build_connectivities())

    return info
