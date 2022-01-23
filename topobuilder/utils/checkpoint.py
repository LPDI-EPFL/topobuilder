# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import Optional, Union, Dict
from pathlib import Path, PosixPath
import json
import sys
from types import GeneratorType

# External Libraries
from logbook import Logger
import numpy as np

# This Library
import topobuilder.core as TBcore

___all__ = ['checkpoint_in', 'checkpoint_out']


class GeneralEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, PosixPath):
            return str(o)
        if isinstance(o, GeneratorType):
            return list(o)
        if isinstance(o, set):
            return list(o)
        if isinstance(o, (np.int_, np.intc, np.intp, np.int8, np.int16, np.int32,
                          np.int64, np.uint8, np.uint16, np.uint32, np.uint64)):
            return int(o)
        elif isinstance(o, (np.float_, np.float16, np.float32, np.float64)):
            return float(o)
        elif isinstance(o, np.ndarray):
            return o.tolist()

        return json.JSONEncoder.default(self, o)


def checkpoint_in( log: Logger, filename: Union[Path, str] ) -> Optional[Dict]:
    """Incoming checkpoint.
    """
    if TBcore.get_option('system', 'forced'):
        return None

    filename = Path(filename)
    if filename.is_file():

        log.info(f'CHECKPOINT: Reloading from {filename}\n')
        with Path(filename).open() as fd:
            try:
                data = json.load(fd)
            except json.JSONDecodeError:
                return None
        return data

    return None


def checkpoint_out( log: Logger, filename: Union[Path, str], data: Dict ):
    """Dump checkpoint.
    """
    filename = Path(filename)

    log.info(f'CHECKPOINT: Creating at {filename}\n')
    with filename.open('w') as fd:
        json.dump(data, fd, cls=GeneralEncoder)
