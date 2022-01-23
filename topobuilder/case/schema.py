# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import getpass
import re
import copy
from collections import OrderedDict

# External Libraries
import numpy as np
from marshmallow import Schema, ValidationError, fields, validates_schema
from marshmallow.validate import Regexp

# This Library
from topobuilder._version import get_versions

__all__ = ['CaseSchema', 'CaseError']

# Default Values
_DEFAULT_HELIX_LENGTH_        = 13
_DEFAULT_BETA_LENGTH_         = 7
_DEFAULT_HELIX_PERIODE_       = 1.5
_DEFAULT_BETA_PERIODE_        = 3.2
_DEFAULT_HELIX_DISTANCE_      = 10
_DEFAULT_HELIX_BETA_DISTANCE  = 11
_DEFAULT_BETA_PAIR_DISTANCE_  = 4.85
_DEFAULT_BETA_STACK_DISTANCE_ = 8
_DEFAULT_LOOP_DISTANCE_       = 18.97
_DEFAULT_LOOP_PERIODE_        = 3.2

_ACCEPTED_SSE_TYPES_          = r'^[HE]$|^S[2-9]\d*$'
_ACCEPTED_SSE_PATTERN_        = re.compile(_ACCEPTED_SSE_TYPES_)
_ACCEPTED_SSE_ERROR_          = "Structure type should meet " \
                                "the pattern: '{}'".format(_ACCEPTED_SSE_TYPES_)

_ACCEPTED_SSE_ID_             = r'^\w\d+[HE]$'
_ACCEPTED_SSE_ID_PATTERN_     = re.compile(_ACCEPTED_SSE_ID_)
_ACCEPTED_SSE_ID_ERROR_       = "Secondary structure id should meet " \
                                "the pattern: '{}'".format(_ACCEPTED_SSE_ID_)


# Global Configuration
class LengthsSchema( Schema ):
    class Meta:
        ordered = True

    H = fields.Integer(default=_DEFAULT_HELIX_LENGTH_,
                       metadata='Number of amino acids in unspecified alpha helix.')
    E = fields.Integer(default=_DEFAULT_BETA_LENGTH_,
                       metadata='Number of amino acids in unspecified beta strand.')


class DistanceSchema( Schema ):
    class Meta:
        ordered = True

    aa = fields.Number(default=_DEFAULT_HELIX_DISTANCE_)
    ab = fields.Number(default=_DEFAULT_HELIX_BETA_DISTANCE)
    bb_pair = fields.Number(default=_DEFAULT_BETA_PAIR_DISTANCE_)
    bb_stack = fields.Number(default=_DEFAULT_BETA_STACK_DISTANCE_)
    max_loop = fields.Number(default=_DEFAULT_LOOP_DISTANCE_)
    loop_step = fields.Number(default=_DEFAULT_LOOP_PERIODE_)

    def get_x_distance( self, data: dict, type1: str, type2: str ) -> float:
        """Provide the x distance between 2 secondary structures depending on their type.
        """
        if type1 is None:
            return 0
        if type1 == 'H' and type2 == 'H':
            return data['aa']
        if type1 == 'H' and type2 == 'E':
            return data['ab']
        if type1 == 'E' and type2 == 'H':
            return data['ab']
        if type1 == 'E' and type2 == 'E':
            return data['bb_pair']

    def get_z_distance( self, data: dict, type1: str, type2: str ) -> float:
        """Provide the z distance between 2 secondary structures depending on their type.
        """
        if type1 is None:
            return 0
        if type1 == 'H' and type2 == 'H':
            return data['aa']
        if type1 == 'H' and type2 == 'E':
            return data['ab']
        if type1 == 'E' and type2 == 'H':
            return data['ab']
        if type1 == 'E' and type2 == 'E':
            return data['bb_stack']


class DefaultSchema( Schema ):
    class Meta:
        ordered = True

    length = fields.Nested(LengthsSchema(), default=LengthsSchema().dump({}))
    distance = fields.Nested(DistanceSchema(), default=DistanceSchema().dump({}))


class ConfigurationSchema( Schema ):
    class Meta:
        ordered = True

    name = fields.String(required=True, default='<name>',
                         error_messages={'required': 'A case identifier is required'},
                         metadata='Case identifier.')
    user = fields.String(default=getpass.getuser(),
                         metadata='User identifier.')
    defaults = fields.Nested(DefaultSchema(), default=DefaultSchema().dump({}),
                             metadata='Default parameters.')
    relative = fields.Boolean(default=True,
                              metadata='Relative vs. absolute coordinates.')
    reoriented = fields.Boolean(default=False,
                                metadata='Has connectivity directions been applied?')
    flip_first = fields.Boolean(default=False,
                                metadata='When applying topologies, start by flipping the first SSE.')
    protocols = fields.List(fields.Dict, default={}, metadata='Pipeline of protocols to run.')
    comments = fields.List(fields.String, default=[get_versions()['version']],
                           metadata='Relative vs. absolute coordinates.')


# Topology
class CoordinateSchema( Schema ):
    class Meta:
        ordered = True

    x = fields.Number(metadata='Value for the X axis.')
    y = fields.Number(metadata='Value for the Y axis.')
    z = fields.Number(metadata='Value for the Z axis.')

    def check_completeness( self, data: dict, info: dict ):
        """Provided de :class:`.CaseSchema` definition is ``absolute`` and not ``relative``,
        :class:`.CoordinateSchema` should have explicitely defined ``x``, ``y`` and ``z``.
        """
        if 'x' not in data:
            raise CaseError('X value for {} of the secondary structures must be provided in absolute mode.'.format(info), 'x')
        if 'y' not in data:
            raise CaseError('Y value for {} of the secondary structures must be provided in absolute mode.'.format(info), 'y')
        if 'z' not in data:
            raise CaseError('Z value for {} of the secondary structures must be provided in absolute mode.'.format(info), 'z')

    def fill_missing( self, data: dict, value: float ) -> dict:
        """Fill non-specified coordinates with the provided value
        """
        for i in ['x', 'y', 'z']:
            data.setdefault(i, value)
        return data

    def append_values( self, data: dict, value: dict ) -> dict:
        """Fill non-specified coordinates with the provided value
        """
        for i in ['x', 'y', 'z']:
            if i in value:
                data[i] += value[i]
        return data

    def distance( self, data1: dict, data2: dict ) -> float:
        """Provide euclidean distance between two coordinates
        """
        return np.linalg.norm(np.asarray([data1['x'], data1['y'], data1['z']]) - np.asarray([data2['x'], data2['y'], data2['z']]))


class StructureSchema( Schema ):
    class Meta:
        ordered = True

    id = fields.String(validate=Regexp(_ACCEPTED_SSE_ID_, error=_ACCEPTED_SSE_ID_ERROR_),
                       metadata='Secondary structure identifier.')
    type = fields.String(required=True, default='<type>',
                         validate=Regexp(_ACCEPTED_SSE_PATTERN_, error=_ACCEPTED_SSE_ERROR_),
                         metadata='Type of secondary structure.')
    length = fields.Integer(metadata='Amino acid length of the secondary structure.')
    coordinates = fields.Nested(CoordinateSchema())
    tilt = fields.Nested(CoordinateSchema())
    layer_tilt = fields.Nested(CoordinateSchema())
    metadata = fields.Dict(metadata='SSE-specific content can be added here by the different plugins.')

    def get_position( self, data: dict ) -> dict:
        """Shortcut to the x, y, z coordinates of the :class:`.StructureSchema`.
        """
        return copy.deepcopy({'x': data['coordinates']['x'], 'y': data['coordinates']['y'], 'z': data['coordinates']['z']})

    def check_completeness( self, data: dict ) -> OrderedDict:
        """Provided de :class:`.CaseSchema` definition is ``absolute`` and not ``relative``,
        :class:`.StructureSchema` should have explicitely defined ``coordinates``, ``angles`` and
        ``length``.
        """
        schema = CoordinateSchema()

        if 'length' not in data:
            raise CaseError('The length of the secondary structures must be provided in absolute mode.', 'length')
        if 'coordinates' not in data:
            raise CaseError('Coordinates of the secondary structures must be provided in absolute mode.', 'coordinates')
        else:
            schema.check_completeness(data['coordinates'], 'coordinates')
        if 'tilt' not in data:
            raise CaseError('Tilt of the secondary structures must be provided in absolute mode.', 'tilt')
        else:
            schema.check_completeness(data['tilt'], 'tilt')

    def cast_absolute( self, data: dict, position: dict, defaults: dict ) -> dict:
        """Transform a ``relative`` :class:`.StructureSchema` into an ``absolute`` one.
        """
        cschema = CoordinateSchema()

        # length
        if 'length' not in data:
            data.setdefault('length', defaults['length'][data['type']])

        # coordinates
        if 'coordinates' not in data:
            data.setdefault('coordinates', {'x': 0, 'y': 0, 'z': 0})
        data['coordinates'] = cschema.fill_missing(data['coordinates'], 0)
        data['coordinates'] = cschema.append_values(data['coordinates'], position)

        # tilt
        if 'tilt' not in data:
            data.setdefault('tilt', {'x': 0, 'y': 0, 'z': 0})
        else:
            data['tilt'] = cschema.fill_missing(data['tilt'], 0)

        # layer_tilt
        if 'layer_tilt' not in data:
            data.setdefault('layer_tilt', {'x': 0, 'y': 0, 'z': 0})
        else:
            data['layer_tilt'] = cschema.fill_missing(data['layer_tilt'], 0)

        # metadata
        data.setdefault('metadata', {})

        return data


class TopologySchema( Schema ):
    class Meta:
        ordered = True

    architecture = fields.List(fields.List(fields.Nested(StructureSchema())),
                               default=[[]], required=True,
                               metadata='Positions and definitions for the secondary structures.')
    connectivity = fields.List(fields.List(fields.String()),
                               metadata='Sequence order of the secondary structures.')


# Case
class CaseSchema( Schema ):
    class Meta:
        ordered = True

    configuration = fields.Nested(ConfigurationSchema(), required=True,
                                  default=ConfigurationSchema().dump({}),
                                  error_messages={'required': 'Configuration data is required'},
                                  metadata='TopoBuilder case definition.')
    topology = fields.Nested(TopologySchema(), required=True, default=TopologySchema().dump({}),
                             error_messages={'required': 'A topological definition is required'},
                             metadata='Topology Definition.')
    metadata = fields.Dict(metadata='Content can be added here by the different plugins.')

    @validates_schema
    def validates_absolute( self, data: dict ):
        """Provided de :class:`.CaseSchema` definition is ``absolute`` and not ``relative``, all
        secondary structures from the defined ``architecture`` should have explicitely defined
        ``coordinates``, ``angles`` and ``length``.
        """
        if data['configuration']['relative'] == False:
            schema = StructureSchema()
            for layer in data['topology']['architecture']:
                for sse in layer:
                    schema.check_completeness(sse)


# Error
class CaseError( ValidationError ):
    """Errors referring to :class:`.CaseSchema` processing"""
