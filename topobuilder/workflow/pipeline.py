# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from copy import deepcopy
from typing import Union, Dict, List
from pathlib import Path
from inspect import signature
import json

# External Libraries
import yaml
import logbook

# This Library
from .log import logger_group, logger_level, LoggedError
from .node import Node, NodeOptionsError, NodeMissingError
from .plugins import TBplugins

__all__ = ['Pipeline', 'ProtocolIncompatibilityError', 'EmptyProtocolError']


class Pipeline( object ):
    """Consecutively executes the :class:`.Nodes` of a protocol.
    """

    def __init__( self, protocol: Union[List[Dict], str, Path] ):
        self.protocol = protocol
        self.nodes = []
        self.current = 0
        self.checked = False
        self.log = logbook.Logger('PIPELINE')
        logger_group.add_logger(self.log)
        logger_level(logger_group)

        if not isinstance(self.protocol, list):
            protocol = str(Path(protocol).resolve())
            try:
                self.protocol = json.loads("".join([x.strip() for x in open(protocol).readlines()]))
                self.log.info(f'Loading JSON protocol from file {protocol}.')
            except json.JSONDecodeError:
                self.protocol = yaml.load(open(protocol), Loader=yaml.Loader)
                self.log.info(f'Loading YAML protocol from file {protocol}.')
        else:
            self.log.info('Loading pre-parsed protocol.')

        # Load the requested plugins for this pipeline
        self.log.debug(f'Protocol contains a total of {len(self.protocol)} nodes.')
        self.log.info('Starting protocol initialization.\n\n')
        for i, node in enumerate(self.protocol):
            self.log.debug(node)
            # Plugin name not specified
            if 'name' not in node:
                err = f'All nodes require a "name" field. None is given for node {i + 1}'
                raise NodeOptionsError(err)
            # Cannot find plugin
            if node['name'] not in TBplugins.source.list_plugins():
                raise NodeMissingError(f'Requested plugin {node["name"]} cannot be found.')
            # Load the plugin class
            this_node = getattr(TBplugins.source.load_plugin(node['name']), node['name'], None)
            # Check it exists
            if this_node is None:
                raise NodeMissingError(f'Requested plugin {node["name"]} does not have a {node["name"]} class.')
            # Check that nodes come from the Node class
            if not issubclass(this_node, Node):
                raise NodeMissingError(f'Class {node["name"]} has to derive from Node.')
            # Load the nodes
            try:
                node = deepcopy(node)
                nname = node.pop('name', None)
                node.pop('tag', None)
                self.nodes.append(this_node(**node, tag=i + 1))
            except TypeError as e:
                err = str(e).replace('__init__()', nname) + '\n'
                err += ' ' * 18 + 'Avaiable keywords are: ' + ', '.join(signature(this_node).parameters.keys())
                raise NodeOptionsError(err)
        self.log.info('-' * 50)

    def check( self, case: Union[Dict, List[Dict]] ) -> List[Dict]:
        """
        """
        self.log.info('Starting protocol checking.\n\n')
        k = [case, ] if not isinstance(case, list) else case
        for node in self.nodes:
            k = node.check(k)
        self.checked = True
        self.log.info('-' * 50)
        return self

    def execute( self, case: Union[Dict, List[Dict]] ) -> List[Dict]:
        """
        """
        self.log.info('Starting protocol execution.\n\n')
        if not self.checked:
            raise UncheckedPipelineError('The pipeline has not been checked before starting running')

        k = [case, ] if not isinstance(case, list) else case
        for node in self.nodes:
            k = node.execute(k)
        self.log.info('-' * 50)
        return k


class UncheckedPipelineError(LoggedError):
    """Raised when trying to run the pipeline without previously check it.
    """


class ProtocolIncompatibilityError(LoggedError):
    """Raised when multiple protocol sources are being called simultaneously.
    """


class EmptyProtocolError(LoggedError):
    """Raised when no protocol data is provided.
    """
