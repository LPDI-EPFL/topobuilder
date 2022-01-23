# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import abc
import sys
import json
from pathlib import PosixPath, Path
from types import GeneratorType
from typing import List, Dict, Union, Optional

# External Libraries
import numpy as np
from logbook import Logger, StreamHandler

# This Library
from .log import logger_group, logger_level, LoggedError

StreamHandler(sys.stdout).push_application()

__all__ = ['Node', 'NodeDataError', 'NodeOptionsError', 'NodeMissingError']


class Node( abc.ABC ):
    """Defines each individual step of the :class:`.Pipeline`.

    This is the base class from which all plugins must derive.
    """
    REQUIRED_FIELDS = ()
    RETURNED_FIELDS = ()
    VERSION = 'v1.0'

    def __init__( self, tag: int ):
        super(Node, self).__init__()
        self.tag = f'{tag:02d}'
        self.nodeID = f'{self.tag}-{type(self).__name__}-{self.VERSION}'
        self.log = Logger(f'Node => {self.nodeID}')
        logger_group.add_logger(self.log)
        logger_level(logger_group)
        self.checked = False
        self.log.info(f'Starting new work node.')

    def check( self, dummy: List[Dict] ) -> List[Dict]:
        """Evaluates the feasability of executing the :class:`.Node`.

        :param dummy: Shape of the data provided by the previous :class:`.Node`.

        :return: Shape of the data after being modified by the current :class:`.Node`.

        :raises:
            :NodeDataError: If the input data shape does not match the expected.
        """
        self.log.info('Checking node viability according to input data.')
        back = []
        self.log.debug(f'Checking a total of {len(dummy):04d} cases.')
        for dt in dummy:
            back.append(self.single_check(dt))
        self.checked = True
        return back

    @abc.abstractmethod
    def single_check( self, dummy: Dict ) -> Dict:
        """Evaluates the feasability of executing the :class:`.Node`.

        :param dummy: Shape of the data provided by the previous :class:`.Node`.

        :return: Shape of the data after being modified by the current :class:`.Node`.

        :raises:
            :NodeDataError: If the input data shape does not match the expected.
        """
        raise NotImplementedError()

    def execute( self, data: List[Dict] ) -> List[Dict]:
        """Process all the data.

        :param data: Data to execute.

        :return: Modified data.
        """
        self.log.info('Executing node.')
        if not self.checked:
            raise NodeUncheckedError(f'{self.nodeID}: Trying to execute when not checked!')
        self.log.debug(f'Executing a total of {len(data):04d} cases.')
        back = []
        for dt in data:
            back.append(self.single_execute(dt))
        self.log.debug(f'Generated a total of {len(back):04d} cases.')
        return back

    @abc.abstractmethod
    def single_execute( self, data: Dict ) -> Dict:
        """Individually process each data entry.

        :param data: Data to execute.

        :return: Modified data.
        """
        raise NotImplementedError()

    def checkpoint_in( self, filename: Union[Path, str] ) -> Optional[Dict]:
        """If a checkpoint exists for this node, load it.

        :param filename: Name of the checkpoint file.

        :return: The recovered data.
        """
        from topobuilder.core import core
        if core.get_option('system', 'forced'):
            return None

        filename = Path(filename)
        if filename.is_file():
            self.log.info(f'CHECKPOINT: Reloading from {filename}')
            with Path(filename).open() as fd:
                try:
                    data = json.load(fd)
                except json.JSONDecodeError:
                    return None
            return data

        return None

    def checkpoint_out( self, filename: Union[Path, str], data: Dict ) -> None:
        """Dump a checkpoint file with the working data.

        :param filename: Name of the checkpoint file.
        :param data: Data to dump.
        """
        filename = Path(filename)
        self.log.info(f'CHECKPOINT: Creating at {filename}')
        with filename.open('w') as fd:
            json.dump(data, fd, cls=GeneralEncoder)


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


class NodeDataError(LoggedError):
    """Raises when the :class:`.Node` is not being provided the rigth data shape.
    """


class NodeOptionsError(LoggedError):
    """Raises when the :class:`.Node` is not provided the necessary options to be executed.
    """


class NodeMissingError(LoggedError):
    """Raises when a requested :class:`.Node` cannot be found or is not the rigth type.
    """


class NodeUncheckedError(LoggedError):
    """Raises when an unchecked :class:`.Node` is requested to be executed.
    """
