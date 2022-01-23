# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import os
import re
import json
import string
import math
import textwrap
from copy import deepcopy
from pathlib import Path
from typing import Optional, Tuple, Dict, List, Union, TypeVar
from collections import OrderedDict
from itertools import chain, zip_longest
from string import ascii_uppercase
import multiprocessing as mp

# External Libraries
import numpy as np
import yaml
try:
    from yaml import CDumper as Dumper
except ImportError:
    from yaml import Dumper
from yaml.representer import SafeRepresenter

# This Library
from .schema import (CaseSchema, CaseError, TopologySchema, CoordinateSchema,
                     StructureSchema, DistanceSchema, ConfigurationSchema)

__all__ = ['Case']

C = TypeVar('C', bound='Case')


class Case( object ):
    """Builds the base Case object.
    """
    def __init__( self, init: Optional[Union[str, dict, Path, C]] = None ):
        self.data = OrderedDict()
        self.schema = CaseSchema()
        if isinstance(init, str):
            self.data = OrderedDict({'configuration': {'name': init}})
        elif isinstance(init, Case):
            self.data = OrderedDict(deepcopy(init.check().data))
        elif isinstance(init, (dict, OrderedDict)):
            self.data = OrderedDict(deepcopy(init))
            self.data = self.schema.load(self.data)
        elif isinstance(init, Path):
            if not init.is_file():
                raise IOError('Unable to find case file {}'.format(init.resolve()))
            if init.suffix == '.gz':
                raise IOError('Unable to manage gzipped file case {}'.format(init.resolve()))
            try:
                self.data = json.loads("".join([x.strip() for x in open(init).readlines()]))
            except json.JSONDecodeError:
                self.data = yaml.load(open(init), Loader=yaml.Loader)

        self.check()

    @property
    def name( self ) -> str:
        """Returns the :class:`.Case` identifier.

        :return: class:`str`
        """
        return self['configuration.name']

    @property
    def shape( self ) -> Tuple[int]:
        """Returns a :class:`tuple` with the number of secondary structures in each layer.

        :return: :class:`tuple`
        """
        if 'architecture' not in self:
            return tuple()
        return tuple([len(layer) for layer in self['topology.architecture']])

    @property
    def shape_len( self ) -> Tuple[Tuple[int]]:
        """Returns the length of each secondary structure in a :class:`tuple` with
        the same shape as the ``architecture``.

        :return: :class:`tuple`
        """
        if 'architecture' not in self:
            return tuple()

        # make a copy absolute first
        c = Case(self.data).cast_absolute()
        out = []
        for layer in c['topology.architecture']:
            out.append([])
            for sse in layer:
                out[-1].append(sse['length'])
            out[-1] = tuple(out[-1])
        return tuple(out)

    @property
    def connectivity_len( self ) -> Tuple[Tuple[int]]:
        """Returns the lengths of the secondary structures in the order of each connectivity.

        :return: :class:`tuple`
        """
        if 'connectivity' not in self:
            return [[]]
        else:
            c = self.cast_absolute()
            lengths = []
            for tplg in c['topology.connectivity']:
                thislen = []
                for sse in tplg:
                    thislen.append(c.get_sse_by_id(sse)['length'])
                lengths.append(thislen)
            return lengths

    @property
    def center_shape( self ) -> Dict:
        """Returns the limit positions for each layer and in total, taking
        only into account the center of each sse.

        :return: :class:`dict`
        """
        if 'architecture' not in self:
            return {}

        asciiU = string.ascii_uppercase

        # make a copy absolute first
        c = Case(self.data).cast_absolute()
        result = {}
        for il, layer in enumerate(c['topology.architecture']):
            result.setdefault(asciiU[il], {'top': 0, 'bottom': 0, 'left': 0, 'right': 0})
            for sse in layer:
                if sse['coordinates']['y'] > result[asciiU[il]]['top']:
                    result[asciiU[il]]['top'] = sse['coordinates']['y']
                if sse['coordinates']['y'] < result[asciiU[il]]['bottom']:
                    result[asciiU[il]]['bottom'] = sse['coordinates']['y']
                if sse['coordinates']['x'] < result[asciiU[il]]['left']:
                    result[asciiU[il]]['left'] = sse['coordinates']['x']
                if sse['coordinates']['x'] > result[asciiU[il]]['right']:
                    result[asciiU[il]]['right'] = sse['coordinates']['x']
            result[asciiU[il]]['width'] = result[asciiU[il]]['right'] - result[asciiU[il]]['left']
            result[asciiU[il]]['hight'] = result[asciiU[il]]['top'] - result[asciiU[il]]['bottom']
        return result

    @property
    def connectivity_count( self ) -> int:
        """Returns the number of connectivities defined in the ``topology``.

        :return: :class:`int`
        """
        if 'connectivity' not in self:
            return 0
        else:
            return len(self['topology.connectivity'])

    @property
    def architecture_str( self ) -> str:
        """Returst a string representation of the architecture.

        :return: :class:`str`
        """
        if 'architecture' not in self:
            return ''
        return architecture_cast(self['topology.architecture'])

    @property
    def connectivities_str( self ) -> Tuple[str]:
        """Returns a list of string representations of topologies.

        :return: :class:`tuple` of :class:`str`
        """
        result = []
        if self.connectivity_count == 0:
            return result

        for i, _ in enumerate(self['topology.connectivity']):
            result.append(topology_cast(self.data, i))
        return tuple(result)

    @property
    def directionality_profile( self ) -> str:
        """Returns a binary-type signature representing the directionality
        of all secondary structures without taking connectivity into account.

        .. note::
            Proper values for this function requires the :class:`Case` to be
            **reoriented**.

        :return: :class:`str`
        """
        result = ''
        if 'architecture' not in self:
            return result

        c = Case(self).cast_absolute()
        for layer in c['topology.architecture']:
            for sse in layer:
                result += '1' if sse['tilt']['x'] > 90 and sse['tilt']['x'] < 270 else '0'
        return result

    @property
    def main_path( self ) -> Path:
        """Main :class:`Path` for the current project.

        :return: :class:`Path`
        """
        return Path(os.path.join(*self.name.split('_')))

    @property
    def undirected_path( self ) -> Path:
        """:class:`Path` to the undirected sketch.

        :return: :class:`Path`
        """
        return Path(os.path.join(self.main_path, 'architecture'))

    @property
    def connectivities_paths( self ) -> Tuple[Path]:
        """Returns a list with the expected :class:`Path` for each connectivity.

        :return: :class:`tuple` of :class:`Path`
        """
        return tuple([Path(os.path.join(self.main_path, 'connectivity', x)) for x in self.connectivities_str])

    @property
    def ordered_structures( self ) -> List[Dict]:
        """Returns the secondary structures in order of connectivity.

        :return: :func:`list` of :class:`dict`
        """
        c = self.cast_absolute()

        if c.is_reoriented:
            if c.connectivity_count != 1:
                raise CaseLogicError('Ordered structures can only be obtained from single-connectivity cases.')
            sse = []
            for s in c.connectivities_str[0].split('.'):
                sse.append(c.get_sse_by_id(s))
            return sse
        else:
            return list(chain(*c['topology.architecture']))

    @property
    def secondary_structure( self ) -> str:
        """DSSP-like secondary structure definition.

        .. note::
            This function requires the :class:`Case` to contain a single connectivity.

            This function requires the existence of ``metadata.loop_lengths``.

        :return: DSSP-like secondary structure definition.

        :raises:
            :CaseLogicError: if ``connectivity_count > 1``.
            :CaseLogicError: if ``metadata.loop_lengths is None``.
            :CaseLogicError: if ``abs(len(sse_list) - len(loop_length)) != 1``.
        """
        c = self.cast_absolute()

        if c.connectivity_count != 1:
            raise CaseLogicError('DSSP string can only be obtained from single-connectivity cases.')

        c = c.apply_topologies()[0]

        loop_length = c['metadata.loop_lengths']
        if loop_length is None:
            raise CaseLogicError('There is no loop length available.')
        loop_length = [''.join(['L', ] * x) for x in loop_length]

        sse_list = [''.join([x['type'], ] * x['length']) for x in c.ordered_structures]
        if abs(len(sse_list) - len(loop_length)) != 1:
            raise CaseLogicError('Number of loops and sse must differ by one.')

        pieces = [sse_list, loop_length] if len(loop_length) < len(sse_list) else [loop_length, sse_list]
        return ''.join([x for x in chain(*zip_longest(*pieces)) if x is not None])

    @property
    def sse_pairing( self ) -> List[str]:
        """Generate Rosetta-like pairing definitions.

        Check Rosetta's **SetSecStructEnergiesMover** for more info.

        .. note::
            This function requires the :class:`Case` to contain a single connectivity.

        :return: Three :class:`str`: ``ss_pair``, ``hh_pair`` and ``hss_triplets``.

        :raises:
            :CaseLogicError: if ``connectivity_count > 1``.
        """
        from topobuilder.case.schema import CoordinateSchema
        c = self.cast_absolute()
        schema = CoordinateSchema()

        if c.connectivity_count != 1:
            raise CaseLogicError('DSSP string can only be obtained from single-connectivity cases.')

        c = c.apply_topologies()[0]

        pfl = c.directionality_profile
        ppfl = {}
        for i, sse in enumerate(c):
            ppfl.setdefault(sse[2]['id'], int(pfl[i]))

        ordsse = c.ordered_structures
        E = [(i + 1, x['id'], ppfl[x['id']]) for i, x in enumerate([x for x in ordsse if x['type'] == 'E'])]
        H = [(i + 1, x['id'], ppfl[x['id']]) for i, x in enumerate([x for x in ordsse if x['type'] == 'H'])]

        ss_pair = []
        SS_pair = []
        for i in range(len(E)):
            for j in range(i + 1, len(E)):
                n1, n2 = E[i][1], E[j][1]
                if n1[0] == n2[0] and abs(int(n1[1]) - int(n2[1])) == 1:
                    ss_pair.append('{}-{}.{}.99'.format(E[i][0], E[j][0], 'A' if E[i][-1] != E[j][-1] else 'P'))
                    SS_pair.append((E[i], E[j]))

        hh_pair = []
        for i, h in enumerate(H):
            for j in range(i + 1, len(H)):
                n1, n2 = h[1], H[j][1]
                ld = abs(ascii_uppercase.find(n1[0]) - ascii_uppercase.find(n2[0]))
                cd = schema.distance(c.get_sse_by_id(n1)['coordinates'], c.get_sse_by_id(n2)['coordinates'])
                if ld == 1 and cd < 15:  # Distance defined in Rosetta's HelixPairingFilter
                        hh_pair.append('{}-{}.{}'.format(h[0], H[j][0], 'A' if h[-1] != H[j][-1] else 'P'))

        hss_triplets = []
        for ss in SS_pair:
            for h in H:
                ld = abs(ascii_uppercase.find(h[1][0]) - ascii_uppercase.find(ss[0][1][0]))
                if ld <= 1:
                    cd1 = schema.distance(c.get_sse_by_id(h[1])['coordinates'],
                                          c.get_sse_by_id(ss[0][1])['coordinates'])
                    cd2 = schema.distance(c.get_sse_by_id(h[1])['coordinates'],
                                          c.get_sse_by_id(ss[1][1])['coordinates'])
                    if cd1 >= 7.5 and cd1 <= 13 and cd2 >= 7.5 and cd2 <= 13:
                        hss_triplets.append('{},{}-{}'.format(h[0], ss[0][0], ss[1][0]))

        return ';'.join(ss_pair), ';'.join(hh_pair), ';'.join(hss_triplets)

    @property
    def potential_connectivities( self ):
        return math.factorial(sum(1 for x in self.architecture_str.split('.')))

    @property
    def potential_sse_orientations( self ):
        from scipy.special import binom
        i = len(self.architecture_str.split('.'))
        return int(binom(i, math.floor(i / 2)) * (1 + i % 2))

    @property
    def is_empty( self ) -> bool:
        return len(self.shape) == 0

    @property
    def is_absolute( self ) -> bool:
        cr = self['configuration.relative']
        if cr is not None:
            return not cr
        return False

    @property
    def is_relative( self ) -> bool:
        return not self.is_absolute

    @property
    def is_reoriented( self ) -> bool:
        cr = self['configuration.reoriented']
        if cr is not None:
            return cr
        return False

    @property
    def flip_first( self ) -> bool:
        """Is the first SSE going to be flipped on topology application?

        :return: :class:`bool`
        """
        return self.cast_absolute()['configuration.flip_first']

    @property
    def mirror_beta_twist( self ) -> bool:
        """Do we need to mirror the corrections to get the right twists?

        :return: :class:`bool`
        """
        return self.cast_absolute()['configuration.mirror_beta_twist']

    def switch_flip_first_to( self, value ) -> C:
        """Change the rule on which SSE to start flipping when applying a topology.

        :param bool value: Values to assign to whether or not to flip the first SSE.

        :return: :class:`.Case`
        """
        c = Case(self)
        c.data['configuration']['flip_first'] = value
        return c

    def get_sse_by_id( self, sse_id: str ) -> Dict:
        """Returns the data corresponding to a given secondary structre according to
        its identifier.

        If the identifier does not belong to any secondary structure of the :class:`.Case`,
        it returns :data:`None`.
        """

        for _, _, sse in self:
            if sse['id'] == sse_id:
                return deepcopy(sse)
        return None

    def get_type_for_layer( self, layer: Union[int, str] ) -> str:
        layerint = layer_int(layer)

        if layerint > len(self['topology.architecture']):
            raise IndexError('Requested layeer is bigger than any available.')

        if len(self['topology.architecture'][layerint]) == 0:
            return 'X'

        return self['topology.architecture'][layerint][0]['type']

    def set_type_for_layer( self, layer: Union[int, str], sse_count: int) -> C:
        sschema = StructureSchema()
        sse_type = self.get_type_for_layer(layer)
        layerint = layer_int(layer)
        layerstr = layer_str(layer)

        ks = Case(self)
        sse = []
        for i in range(sse_count):
            sse.append(sschema.dump({'id': '{0}{1}{2}'.format(layerstr, i + 1, sse_type),
                                     'type': sse_type}))
        ks.data['topology']['architecture'][layerint] = sse
        return ks

    def check( self ) -> C:
        """Evaluate the :class:`.Case` content thourhg the :class:`.CaseSchema`.
        """
        self.data = self.schema.load(self.schema.dump(self.data))
        return self

    def add_architecture( self, architecture: Optional[str] = None ) -> C:
        """Adds an architecture definition to the :class:`.Case`.

        According to `CATH <http://www.cathdb.info/wiki/doku/?id=faq#what_is_cath>`_, an
        architecture defines *structures that are classified according to their overall
        shape as determined by the orientations of the secondary structures in 3D space
        but ignores the connectivity between them*.

        In the **TopoSuite**, we have adopted and adapted this definition to a layer-based
        `FORM <https://doi.org/10.1016/j.str.2009.07.012>`_ definition of secondary structure
        placements.

        The format of an architecture string definition is as follows (lower or upper case
        are both accepted)::

            2h.4e.2h

        This defines a 3-layer topology (each layer separated by points) formed by a first
        layer with 2 helices, a mid-layer with 4 beta strands and a final layer with 2 helices.

        Normally, the length of each secondary structure will be determined by the appropriate
        default setting. If length want to be defined, this can be done by providing length data
        to each layer using `:`::

            2h:13:10.4e:5:5:5:5.2h:7:13

        Notice that, if secondary structure residue length is provided, **it has to be provided
        for all secondary structures**. Unexpected behaviour might arise from doing otherwise.

        .. note::
            All architecture string definitions can be processed by the **TopoSuite**, but not
            all constructs generated by the **TopoSuite** can be minimized into an architecture
            string definition.

        :param str architecture: Architecture string definition.

        :return: :class:`.Case` - updated case.

        :raises:
            :CaseOverwriteError: If an architecture is already defined in this :class:`.Case`
            :CaseError: If the string cannot be properly parsed.

        .. seealso::
            :func:`.add_topology`
        """
        if architecture is None:
            return Case(self.data).check()

        if 'architecture' in self:
            raise CaseOverwriteError('An arquitecture is already defined.')

        c = Case(self.data)
        c.data['topology']['architecture'] = architecture_cast(architecture)['architecture']
        return c.check()

    def add_topology( self, topology: Optional[str] = None ) -> C:
        """Adds a topology definition to the :class:`.Case`.

        According to `CATH <http://www.cathdb.info/wiki/doku/?id=faq#what_is_cath>`_, a
        topology defines *structures are grouped into fold groups at this level depending
        on both the overall shape and connectivity of the secondary structures*.

        In the **TopoSuite**, a topology string defines the connectivity of the secondary
        structures and, by using the `FORM <https://doi.org/10.1016/j.str.2009.07.012>`_ systematic
        naming system (a ``row==letters``, ``columns==numbers`` grid system), it also defines their
        position. As such, a definition such as this::

            B2E.C1H.B1E.A1H.B3E.A2H.B4E.C2H

        would represent a 3-layered topology (layers A, B and C) with two structures in the first
        and third layers and 4 in the middle one.

        By default, residue length of each secondary structure is defined by the default configuration,
        but it can be set up **on an individual basis** by providing the length after the secondary
        structure type::

            B2E5.C1H12.B1E4.A1H16.B3E5.A2H13.B4E5.C2H13

        :param str topology: Topology string definition.

        :return: :class:`.Case` - updated case.

        :raises:
            :CaseError: If the string cannot be properly parsed.
            :CaseOverwriteError: If an architecture is already defined in this :class:`.Case` and the
                provided topology does not match

        .. seealso::
            :func:`.describe_architecture`
        """
        if topology is None:
            return Case(self.data).check()

        t = Case('temp')
        t.data['topology'] = topology_cast(topology)

        if 'architecture' in self:
            if self.architecture_str != t.architecture_str or self.shape_len != t.shape_len:
                raise CaseOverwriteError('Provided topology does not match existing architecture.')

        c = Case(self.data)
        if 'connectivity' not in c:
            c.data['topology'] = topology_cast(topology)
        else:
            conn = topology_cast(topology)['connectivity'][0]
            if ".".join(conn) not in set(self.connectivities_str):
                c.data['topology']['connectivity'].append(conn)
        return c.check()

    def add_secured_topologies( self, topologies: List[C] ) -> C:
        """Adds topologies with their respective connectivity.
        """
        c = Case(self.data)
        c.data['topology'].setdefault('connectivity', [])
        c.data['topology']['connectivity'] = topologies
        return c

    def cast_absolute( self ) -> C:
        """Transform a ``relative`` :class:`.CaseSchema` into an ``absolute`` one.
        """
        c = Case(self.data)
        if self.is_absolute:
            return c

        if c.is_empty:
            raise CaseLogicError('An empty case cannot be made absolute.')

        c.data['configuration']['relative'] = False
        sschema = StructureSchema()
        dschema = DistanceSchema()

        position = {'x': 0, 'y': 0, 'z': 0}
        defaults = c['configuration.defaults']
        zlayer = 0
        for i, layer in enumerate(c['topology.architecture']):
            position['x'] = 0
            back = None if i == 0 else c['topology.architecture'][i - 1][0]['type']
            here = layer[0]['type']
            zlayer = dschema.get_z_distance(defaults['distance'], back, here) * i
            for j, sse in enumerate(layer):
                left = None if j == 0 else layer[j - 1]['type']
                here = sse['type']
                # X shift is inherited in the following structure.
                position['x'] += dschema.get_x_distance(defaults['distance'], left, here)
                position['y'] = 0  # Reset Y coordinate
                position['z'] = zlayer  # Reset Z coordinate
                c.data['topology']['architecture'][i][j] = sschema.cast_absolute(sse, position, defaults)
                position = sschema.get_position(c['topology.architecture'][i][j])

        return c.check()

    def apply_topologies( self ) -> List[C]:
        """Generates a :class:`List` of :class:`.Case` in which the different available connectivities
        have been applied.

        Will run async pool when more than 100 connectivities are provided.

        :return: :class:`List` of :class:`.Case`

        :raises:
            :CaseIncompleteError: If there is not enough data to apply topologies.
        """
        if 'connectivity' not in self or 'architecture' not in self:
            raise CaseIncompleteError('Unable to apply corrections to non-existing fields.')

        cn = self['topology.connectivity']

        if len(cn) == 1:
            if self.is_reoriented:
                return [self, ]
            return [make_topology(self, 0), ]

        if len(cn) > 1:
            if self.is_reoriented:
                raise CaseLogicError('The case has multiple connectivities but its labeled as oriented?')
            if len(cn) < 100:
                return [make_topology(self, i) for i in range(len(cn))]
            else:
                result = []
                pool = mp.Pool(3)
                for i in range(len(cn)):
                    result.append(pool.apply_async(make_topology, args=(self, i)))
                pool.close()
                pool.join()
                return [r.get() for r in result]

    def apply_corrections( self, corrections: Optional[Union[Dict, str, Path]] = None ) -> C:
        """Applies the corrections to the :term:`SKETCH`.

        :return: corrected :class:`.Case`
        """
        if corrections is None:
            return Case(self.data).check()

        if isinstance(corrections, str):
            corrections = Path(corrections)

        if isinstance(corrections, Path):
            if not corrections.is_file():
                raise IOError('Unable to find corrections file {}'.format(corrections.resolve()))
            if corrections.suffix == '.gz':
                raise IOError('Unable to manage gzipped file case {}'.format(corrections.resolve()))
            try:
                crr = json.loads("".join([x.strip() for x in open(corrections).readlines()]))
            except json.JSONDecodeError:
                crr = yaml.load(open(corrections), Loader=yaml.Loader)
            return Case(self.data).apply_corrections(crr)

        if isinstance(corrections, dict) and not bool(corrections):
            return Case(self.data).check()

        # APPLY CONFIGURATION CORRECTIONS (before cast absolute is applied in layer_corrections)
        c = Case(self)
        c = configuration_corrections(corrections, c)

        # APPLY LAYER CORRECTIONS
        corrections = layer_corrections(corrections, c)

        # APPLY SSE CORRECTIONS
        return sse_corrections(corrections, c)

    def set_protocol_done( self, protocol_id: int ) -> C:
        """Label a protocol as executed.

        :param int protocol_id: Identifier of the protocol according to its position.
        """
        if protocol_id == -1:
            return self
        if self['configuration.protocols'] is None or len(self['configuration.protocols']) < protocol_id:
            raise IndexError('Trying to access an unspecified protocol.')

        c = Case(self.data)
        c.data['configuration']['protocols'][protocol_id].setdefault('status', True)
        c.data['configuration']['protocols'][protocol_id]['status'] = True
        return c

    def assign_protocols( self, protocols: List[Dict] ) -> C:
        """Overwrite current protocols in the :class:`.Case` with the provided ones.

        :param protocols: New protocols for the :class:`.Case`
        """
        c = Case(self.data)
        if c['configuration.protocols'] is None:
            c.data['configuration'].setdefault('protocols', [])
        c.data['configuration']['protocols'] = protocols
        return c

    def write( self, prefix: Optional[Union[str, Path]] = None, format: str = 'yaml' ) -> Path:
        """Write :class:`.Case` into a file.

        :param str prefix: If not specified, output file is created in the **current working directory**,
            using :class:`.Case` ``configuration.name`` as prefix. A :class:`str` will modify the name of
            the file. A directory :class:`Path` will change dir and use the default prefix, while a
            non-directory :class:`Path` will set up output directory and file prefix.
        :type prefix: Union[str, Path]
        :param str format: Output format. Options are [``yaml``, ``json``]

        :return: :class:`Path` - filename

        :raises:
            :ValueError: If format is not ``yaml`` or ``json``.

        """
        if format not in ['yaml', 'json']:
            raise ValueError('Available formats are yaml or json.')

        if prefix is None:
            prefix = self['configuration.name']
        if isinstance(prefix, Path):
            if prefix.is_dir():
                prefix = prefix.joinpath(self['configuration.name'])

        if format == 'yaml':
            YAML_Dumper()
            with open('{}.yml'.format(str(prefix)), 'w') as fd:
                yaml.dump(self.data, fd, Dumper=Dumper, default_flow_style=False)
            return Path('{}.yml'.format(str(prefix)))
        else:
            with open('{}.json'.format(str(prefix)), 'w') as fd:
                json.dump(self.data, fd, indent=2)
            return Path('{}.json'.format(str(prefix)))

    def __contains__( self, item ):
        if isinstance(item, str):
            if item == 'architecture':
                ta = self['topology.architecture']
                return False if ta is None else ta != [[]]
            if item == 'connectivity':
                ta = self['topology.connectivity']
                return False if ta is None else ta != [[]]

        raise NotImplementedError()

    def __iter__( self ):
        architecture = self['topology.architecture']
        for i, layer in enumerate(architecture):
            for j, sse in enumerate(layer):
                yield i, j, sse

    def __getitem__( self, key ):
        r = self.data
        for k in key.split('.'):
            r = r.get(k, None)
            if r is None:
                break
        return r

    def __eq__( self, value ):
        if not isinstance(value, Case):
            raise NotImplementedError()
        return self.data == value.data

    def __len__( self ):
        return int(np.sum(np.asarray(self.shape)))

    def _ipython_display_( self ):
        from IPython.core.display import display, HTML
        from jinja2 import Template
        sse = {'H': 0, 'E': 0}

        if not self.is_empty:
            for x in self.architecture_str.split('.'):
                sse.setdefault(x[-1], 0)
                sse[x[-1]] += int(x[:-1])
        tplt = Template(textwrap.dedent("""\
            <table class="case_table" style="font-family:monospace;">
                <tr><th style="text-align:center;" colspan="3"><b>{{n}}</b></th></tr>
                <tr>
                    <td style="text-align:left;" colspan="2"><b>Layers</b></td>
                    <td style="text-align:center;">{{y}}</td>
                </tr>
                <tr>
                    <td style="text-align:left;" rowspan="2"><b>SSE</b></td>
                    <td style="text-align:center;"><b>H</b></td>
                    <td style="text-align:center;">{{H}}</td>
                </tr>
                <tr>
                    <td style="text-align:center;"><b>E</b></td>
                    <td style="text-align:center;">{{E}}</td>
                </tr>
                <tr>
                    <td style="text-align:left;" rowspan="2"><b>Folds</b></td>
                    <td style="text-align:center;"><b>Potential</b></td>
                    <td style="text-align:center;">{{f}}</td>
                </tr>
                <tr>
                    <td style="text-align:center;"><b>Specified</b></td>
                    <td style="text-align:center;">{{s}}</td>
                </tr>
            </table>""")).render(n=self.name, **sse, y=len(self.shape),
                                 f=math.factorial(sum(sse[x] for x in sse)),
                                 s=self.connectivity_count)
        display(HTML(tplt))


# Function is global so that it can be pickled -> apply_async
def make_topology( case, count ):
    """
    """

    c = Case(case.data)
    corrections = {}
    c['topology']['connectivity'] = [c['topology.connectivity'][count]]
    for turn in c['topology.connectivity'][0][1 if not c.flip_first else 0::2]:
        corrections.setdefault(turn, {'tilt': {'x': 180}})
    c = c.apply_corrections(corrections)
    c.data['configuration']['reoriented'] = True
    return c


def configuration_corrections( corrections: dict, case: Case) -> Dict:
    """Apply specific configuration-related changes to the :class:`.Case`.

    :param corrections: Existing corrections.
    :param case: Data to guide the corrections.

    :return: Updated corrections.
    """
    cr = corrections.pop('configuration', None)
    if cr is None:
        return case

    cfg = case.data['configuration']
    for k in cr:
        cfg[k] = cr[k]

    schema = ConfigurationSchema()
    case.data['configuration'] = schema.dump(cfg)

    return case


def layer_corrections( corrections: dict, case: Case) -> Dict:
    """Transform layer-specified corrections into SSE ones that can be translated
    into proper :class:`.Case` definitions.

    :param corrections: Existing corrections.
    :param case: Data to guide the corrections.

    :return: Updated corrections.
    """
    case = case.cast_absolute()
    asciiU = string.ascii_uppercase
    sizes = case.center_shape
    maxwidth = max(sizes[l]['width'] for l in sizes)
    for j, layer in enumerate(case['topology.architecture']):
        layer_id = asciiU[j]
        if layer_id in corrections:
            c = corrections[layer_id]
            if 'xalign' in c:
                diffw = maxwidth - sizes[layer_id]['width']
                corrections = lc_xalign(corrections, layer, c['xalign'].lower(), diffw)
            if 'yalign' in c:
                corrections = lc_yalign(corrections, layer, c['yalign'].lower(), case)
            if 'zcurve' in c:
                corrections = lc_zcurve(corrections, layer, c['zcurve'])
    return corrections


def lc_xalign( corrections: Dict, layer: List[Dict], xalign: str, diffw: float ) -> Dict:
    """Apply xalign layer corrections.

    :param corrections: Existing corrections to where the new ones will be accomulated.
    :param layer: List of SSE elements in the layer.
    :param xalign: Value of the ``xalign`` correction. Can be [``left``, ``center``, ``rigth``]
    :param diffw: Difference between the wider layer and the current one.

    :return: Updated corrections.

    This layer correction is applied to the first SSE element of the layer, as the position in the x
    axis of the rest are calculated from this one.
    """
    # Do not apply when (A) this is the widest layer or (B) alignment is left (default)
    if diffw == 0 or xalign == 'left':
        return corrections

    sse_id = layer[0]['id']
    corrections.setdefault(sse_id, {}).setdefault('coordinates', {}).setdefault('x', 0)
    if xalign == 'right':
        corrections[sse_id]['coordinates']['x'] += diffw
    elif xalign == 'center':
        corrections[sse_id]['coordinates']['x'] += diffw / 2
    return corrections


def lc_yalign( corrections: Dict, layer: List[Dict], yalign: str, case: Case ) -> Dict:
    """Apply yalign layer corrections.

    :param corrections: Existing corrections to where the new ones will be accomulated.
    :param layer: List of SSE elements in the layer.
    :param yalign: Value of the ``yalign`` correction. Can be [``top``, ``middle``, ``bottom``]
    :param case: Full data to get the SSE hights.

    :return: Updated corrections.

    This layer correction is applied to each SSE element of the layer. If the SSE lengths are changed
    after this correction is applied, the final content might not be as expected.
    """
    # Do not apply when (A) there is only one SSE or (B) alignemnt is middle (default)
    if len(layer) == 1 or yalign == 'middle':
        return corrections

    tops, bottoms = layer_hights(case, layer)
    for i, sse in enumerate(layer):
        topdiff = abs(max(tops) - tops[i])
        bottomdiff = abs(min(bottoms) - bottoms[i])
        corrections.setdefault(sse['id'], {}).setdefault('coordinates', {}).setdefault('y', 0)
        if yalign == 'top':
            corrections[sse['id']]['coordinates']['y'] += topdiff
        elif yalign == 'bottom':
            corrections[sse['id']]['coordinates']['y'] -= bottomdiff
    return corrections


def lc_zcurve( corrections: Dict, layer: List[Dict], zcurve: float ) -> Dict:
    """
    """
    def define_circle(p1, p2, p3):
        """Returns the center and radius of the circle passing the given 3 points.
        In case the 3 points form a line, returns (None, infinity).
        """
        temp = p2[0] * p2[0] + p2[1] * p2[1]
        bc = (p1[0] * p1[0] + p1[1] * p1[1] - temp) / 2
        cd = (temp - p3[0] * p3[0] - p3[1] * p3[1]) / 2
        det = (p1[0] - p2[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p2[1])

        if abs(det) < 1.0e-6:
            return (None, np.inf)

        # Center of circle
        cx = (bc * (p2[1] - p3[1]) - cd * (p1[1] - p2[1])) / det
        cy = ((p1[0] - p2[0]) * cd - (p2[0] - p3[0]) * bc) / det

        radius = np.sqrt((cx - p1[0])**2 + (cy - p1[1])**2)
        return ((cx, cy), radius)

    def angle_points(a, b, c):
        ba = np.asarray(a) - np.asarray(b)
        bc = np.asarray(c) - np.asarray(b)

        cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
        return np.arccos(cosine_angle)

    def point_on_circle(center, radius, angle):
        """Finding the x, y coordinates on circle, based on given angle
        """
        x = center[0] + (radius * np.cos(angle))
        y = center[1] + (radius * np.sin(angle))

        return x, y

    # Do not apply when (A) there is less than SSE or (B) curve is 0 (default)
    if len(layer) <= 2 or zcurve == 0:
        return corrections

    # Create curve guiding points
    centerINI = layer[0]['coordinates']
    centerINI = [centerINI['x'], centerINI['z']]
    centerEND = layer[-1]['coordinates']
    centerEND = [centerEND['x'], centerEND['z']]
    midEND = (np.asarray(centerEND) - np.asarray(centerINI)) / 2
    centerINI[-1] += zcurve
    centerEND[-1] += zcurve
    midEND[-1] -= zcurve
    points = np.asarray([centerINI, midEND, centerEND])
    center, radius = define_circle(*points)

    N = len(layer)
    point_0 = point_on_circle(center, radius, 0)
    leftmost_angle = -angle_points(points[0], center, point_0)
    rightmost_angle = -angle_points(points[-1], center, point_0)
    step_angle = -(leftmost_angle - rightmost_angle) / (N - 1)

    points = []
    for x in range(N):
        points.append(point_on_circle(center, radius, leftmost_angle + (step_angle * x)))

    for i, sse in enumerate(layer):
        corrections.setdefault(sse['id'], {}).setdefault('coordinates', {}).setdefault('z', 0)
        corrections[sse['id']]['coordinates'].setdefault('x', 0)
        corrections[sse['id']]['coordinates']['z'] += points[i][1]
        if i == 0 or i == len(layer) - 1:
            continue
        corrections[sse['id']]['coordinates']['x'] += points[i][0] - sse['coordinates']['x']

    # # get x and y vectors
    # x = points[:, 0]
    # y = points[:, 1]
    #
    # # calculate polynomial
    # z = np.polyfit(x, y, 2)
    # f = np.poly1d(z)
    #
    # for i, sse in enumerate(layer):
    #     corrections.setdefault(sse['id'], {}).setdefault('coordinates', {}).setdefault('z', 0)
    #     currentX = sse['coordinates']['x']
    #     corrections[sse['id']]['coordinates']['z'] += np.around(f(currentX), decimals=3)

    return corrections

    # # Do not apply when (A) there is less than SSE or (B) curve is 0 (default)
    # if len(layer) <= 2 or zcurve == 0:
    #     return corrections
    #
    # # Create curve guiding points
    # centerINI = layer[0]['coordinates']
    # centerINI = [centerINI['x'], centerINI['z']]
    # centerEND = layer[-1]['coordinates']
    # centerEND = [centerEND['x'], centerEND['z']]
    # midEND = (np.asarray(centerEND) - np.asarray(centerINI)) / 2
    # centerINI[-1] += zcurve
    # centerEND[-1] += zcurve
    # midEND[-1] -= zcurve
    # points = np.asarray([centerINI, midEND, centerEND])
    # # get x and y vectors
    # x = points[:, 0]
    # y = points[:, 1]
    #
    # # calculate polynomial
    # z = np.polyfit(x, y, 2)
    # f = np.poly1d(z)
    #
    # for i, sse in enumerate(layer):
    #     corrections.setdefault(sse['id'], {}).setdefault('coordinates', {}).setdefault('z', 0)
    #     currentX = sse['coordinates']['x']
    #     corrections[sse['id']]['coordinates']['z'] += np.around(f(currentX), decimals=3)
    #
    # return corrections


def sse_corrections( corrections: dict, case: Case ) -> Case:
    """
    """
    cs = CoordinateSchema()
    case = deepcopy(case.data)
    for j, layer in enumerate(case['topology']['architecture']):
        for i, sse in enumerate(layer):
            if sse['id'] in corrections:
                for c in corrections[sse['id']]:
                    if not isinstance(corrections[sse['id']][c], (dict, OrderedDict)):
                        case['topology']['architecture'][j][i].setdefault(c, None)
                        case['topology']['architecture'][j][i][c] = corrections[sse['id']][c]
                    else:  # has to be in ['coordinates', 'tilt', 'layer_tilt']
                        case['topology']['architecture'][j][i].setdefault(c, {})
                        ks = cs.fill_missing(case['topology']['architecture'][j][i][c], 0)
                        case['topology']['architecture'][j][i][c] = cs.append_values(ks, corrections[sse['id']][c])
    return Case(case).check()


def layer_cast( layer: Union[int, str] ) -> Union[int, str]:
    if isinstance(layer, int):
        return string.ascii_uppercase[layer]
    elif isinstance(layer, str):
        return string.ascii_uppercase.find(layer.upper())
    else:
        raise ValueError('Layer is defined by integer or string.')


def layer_int( layer: Union[str, int] ) -> int:
    if isinstance(layer, int):
        return layer
    else:
        return layer_cast(layer)


def layer_str( layer: Union[str, int] ) -> str:
    if isinstance(layer, str):
        return layer.upper()
    else:
        return layer_cast(layer)


def layer_hights( case: Case, layer: List[Dict] ) -> List[List[float]]:
    """
    """
    from .schema import _DEFAULT_BETA_PERIODE_, _DEFAULT_HELIX_PERIODE_

    sc = CoordinateSchema()
    tops = []
    bots = []
    defaults = case['configuration.defaults.length']
    for i, sse in enumerate(layer):
        length = sse['length'] if 'length' in sse else defaults[sse['type']]
        periode = _DEFAULT_BETA_PERIODE_ if sse['type'] == 'E' else _DEFAULT_HELIX_PERIODE_
        max_dist = float(periode * (length - 1))
        centre = sc.fill_missing(sse['coordinates'], 0) if 'coordinates' in sse else sc.fill_missing({}, 0)
        centre = [centre['x'], centre['y'], centre['z']]
        tops.append((np.copy(centre) + np.array([0, max_dist / 2, 0]))[1])
        bots.append((np.copy(centre) - np.array([0, max_dist / 2, 0]))[1])
    return tops, bots


def architecture_cast( architecture: Union[str, List, dict] ) -> Union[dict, str]:
    """
    """
    if isinstance(architecture, str):
        expression = re.compile(r'^(\d+)([EH])$')
        architecture = architecture.upper()
        asciiU = string.ascii_uppercase
        result = {'architecture': []}

        for layer in architecture.split('.'):
            layer = layer.split(':')
            m = re.match(expression, layer[0])
            if not m:
                raise CaseError('Architecture format not recognized.')
            result['architecture'].append([])
            for i in range(int(m.group(1))):
                name = '{0}{1}{2}'.format(asciiU[len(result['architecture']) - 1],
                                          i + 1, m.group(2))
                result['architecture'][-1].append({'type': m.group(2), 'id': name})
                if len(layer) > 1:
                    try:
                        result['architecture'][-1][-1].setdefault('length', int(layer[i + 1]))
                    except IndexError:
                        print('Lengths were not provided for all defined secondary structures.')
                    except ValueError:
                        print('Length values MUST BE integers.')
                    except Exception as e:
                        print(e)

        schema = TopologySchema()
        return schema.dump(result)

    if isinstance(architecture, list):
        architecture = {'architecture': architecture}

    result = []
    for layer in architecture['architecture']:
        result.append('{0}{1}'.format(len(layer), layer[0]['type']))
    result = '.'.join(result)
    return result


def topology_cast( topology: Union[str, dict], count: Optional[int] = 0 ) -> Union[dict, str]:
    """
    """
    if isinstance(topology, str):
        expression = re.compile(r'^([A-Z])(\d+)([EH])(\d*)$')
        topology = topology.upper()
        result = {'architecture': [], 'connectivity': []}
        architecture = []
        for topo in topology.split(','):
            tp = {}
            result['connectivity'].append([])
            architecture.append([])
            for sse in topo.split('.'):
                m = re.match(expression, sse)
                if not m:
                    raise CaseError('Topology format not recognized.')
                sse_id = '{0}{1}{2}'.format(m.group(1), m.group(2), m.group(3))
                result['connectivity'][-1].append(sse_id)
                tp.setdefault(string.ascii_uppercase.find(m.group(1)) + 1,
                              {}).setdefault(int(m.group(2)), (m.group(3), sse_id, m.group(4)))

            if list(sorted(tp.keys())) != list(range(min(tp.keys()), max(tp.keys()) + 1)):
                raise CaseError('Topology format skips layers.')
            for k1 in sorted(tp.keys()):
                architecture[-1].append([])
                if list(sorted(tp[k1].keys())) != list(range(min(tp[k1].keys()), max(tp[k1].keys()) + 1)):
                    raise CaseError('Topology format skips positions in layer {}.'.format(k1))
                for k2 in sorted(tp[k1].keys()):
                    scaffold = {'type': tp[k1][k2][0], 'id': tp[k1][k2][1]}
                    if tp[k1][k2][2] != '':
                        scaffold.setdefault('length', tp[k1][k2][2])
                    architecture[-1][-1].append(scaffold)
        arch_str = []
        for a in architecture:
            astr = []
            for l in a:
                [astr.append(ss['id']) for ss in l]
            arch_str.append(''.join(astr))
        arch_str = list(set(arch_str))
        if len(arch_str) > 1:
            raise CaseLogicError('A case can only contain one architecture.')
        else:
            result['architecture'] = list([list(x) for x in architecture[0]])
        schema = TopologySchema()
        return schema.dump(result)

    return ".".join(topology['topology']['connectivity'][count])


def YAML_Dumper():
    # This is required for YAML to properly print the Schema as an OrderedDict
    # Adapted from https://gist.github.com/oglops/c70fb69eef42d40bed06 to py3
    def dict_representer(dumper, data):
        return dumper.represent_dict(data.items())

    Dumper.add_representer(OrderedDict, dict_representer)
    Dumper.add_representer(str, SafeRepresenter.represent_str)


# Errors
class CaseOverwriteError( CaseError ):
    """Error raised when trying to override :class:`.Case` data in an unexpected way
    """


class CaseIncompleteError( CaseError ):
    """Errors raised when trying to execute a process without enough information
    """


class CaseLogicError( CaseError ):
    """Errors raised when there are incoherent :class:`.Case conditions`
    """
