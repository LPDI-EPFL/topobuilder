# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import sys

# External Libraries
import logbook

# This Library
from topobuilder.core import core

__all__ = ['logger_group', 'logger_level', 'LoggedError']

logbook.StreamHandler(sys.stdout).push_application()

logger_group = logbook.LoggerGroup()


def logger_level( logger_group ):
    if core.get_option('system', 'debug'):
        logger_group.level = logbook.DEBUG
    elif core.get_option('system', 'verbose'):
        logger_group.level = logbook.INFO
    else:
        logger_group.level = logbook.NOTICE


class LoggedError( Exception ):
    """Adds the error message to the log.
    """
    def __init__( self, *args, **kwargs ):
        Exception.__init__(self, *args, **kwargs)
        self.log = logbook.Logger(f'{type(self).__name__}')
        logger_group.add_logger(self.log)
        self.log.error(self)
