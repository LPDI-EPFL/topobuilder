# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>
.. codeauthor:: Zander Harteveld <zandermilanh@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import List, Union, Dict, Optional
from pathlib import Path
from collections import namedtuple
import re

# External Libraries
import pandas as pd

# This Library
from topobuilder.case import Case
from topobuilder.utils import reverse_motif, StructuralError
from topobuilder.case.schema import _ACCEPTED_SSE_ID_
from topobuilder.workflow import Node, NodeOptionsError, NodeDataError


__all__ = ['motif_picker']

#Motif = namedtuple('Motif', ['motifs', 'binder', 'hotpsots', 'attach', 'selection', 'identifier'])


class motif_picker( Node ):
    """Recovers a motif of interest from a protein structure to map upon a :term:`FORM`.

    .. note::
        By adding a motif linking it to one or more secondary structres, several restrictions are added to the system.
        These may include relative positioning, sequence order and positioning, amongst others.

    .. admonition:: To Developers

        A part from adding the ``metadata.motif_picker`` field to the general :class:`.Case`, this :class:`.Node` adds
        information to the ``metadata`` field of selected secondary structures. Ignoring these metadata in a
        :class:`.Node` that alters positioning and orientation of the secondary structures will result in broken
        functional motifs. When applying corrections it should be done through the :class:`.Node` type :class:`.corrector`,
        in which this restrictions are already taken into account.

    :param source: Path to the PDB file of interest containing the motif.
    :param selection: Range selection of the motif segments. Position **must be defined as chain and PDB identifier**,
        that means number and insertion code, if any. Thus, one motif could be selected as ``A:10-25`` when belonging
        to chain ``A``, while multi-segment motifs are made of selections separated by ``,``.
    :param hotspot: Selection of a single position that belongs to the exposed part of the interface. This residue is
        used as a guide to understand the motif's placement. This selection has to follow the same format as selection,
        but without a range.
    :param attach: :term:`FORM` name of the secondary structure to which the segments need to be attached.
        If there is more than one segment, names should be provided in a list.
    :param binder: Range selection for the binder. It can be provided as a range same as ``selection`` or as one or a
        comma-separated list of chain identifiers.
    :param identifier: Name to give the motif. Otherwise a name is picked related to the position of the :class:`.Node`
        in the :class:`.Pipeline`.

    :raises:
        :NodeOptionsError: On **initialization**. If the provided structure source file does not exist.
        :NodeOptionsError: On **initialization**. If the number of requested segments does not match the
            number of secondary structures to attach them.
        :NodeOptionsError: On **initialization**. If the number of requested segments does not match the
            shape definition.
        :NodeOptionsError: On **initialization**. If any secondary structure identifier in ``attach`` does
            not match *TopoBuilder* naming system.
        :NodeDataError: On **check**. If the required fields to be executed are not there.

    """
    SHAPE_PATTERN = r'^[HE]([lx][HE])*$'
    REQUIRED_FIELDS = ('topology.architecture', )
    RETURNED_FIELDS = ('metadata.motif_picker.motifs')
    VERSION = 'v1.0'

    def __init__( self, tag: int,
                  source: Union[Path, str],
                  selection: str,
                  hotspot: str,
                  attach: Union[str, List[str]],
                  binder: Optional[str] = None,
                  identifier: Optional[str] = None ):
        super(motif_picker, self).__init__(tag)

        self.Motif = namedtuple('Motif', ['motifs', 'binder', 'hotpsots', 'attach', 'selection', 'identifier'])

        self.identifier = identifier if identifier is not None else self.tag
        # Obtaining structure input file
        self.source = source.resolve() if isinstance(source, Path) else Path(source).resolve()
        if not self.source.is_file():
            raise NodeOptionsError(f'Provided structure file {self.source} cannot be found.')

        # Check that the number of selection, attach and shape make sense
        self.selection = selection.split(',')
        self.hotspot = hotspot
        self.attach = attach if isinstance(attach, list) else [attach, ]
        #self.shape = shape if isinstance(shape, list) else [shape, ]

        # Each segments belongs to a SSE
        #if len(self.selection) != len(self.attach):
        #    err = 'The number of selection segments should match the number of structures to attach. There are '
        #    err += f'{len(self.selection)} segments requested to be assigned to {len(self.attach)} structures.'
        #    raise NodeOptionsError(err)

        # The names of the SSE to attach the motif can exist.
        err = []
        for sse in self.attach:
            if not re.match(_ACCEPTED_SSE_ID_, sse):
                err.append(f'Unrecognized SSE identifier in {sse}.')
        if len(err) > 0:
            raise NodeOptionsError('\n'.join(err))

        # Capture binder selection
        self.binder = binder

        # Chatty
        self.log.debug(f'Picking a {len(self.selection)}-segment motif from {self.source}.')
        self.log.debug(f'Hotspot residue used to guide orientation is {self.hotspot}.')
        self.log.debug(f'Motif segments are to be assigned to SSE(s) {",".join(self.attach)}.')
        if self.binder is not None:
            self.log.debug(f'A binder on selection {self.binder} has been included.')


    def single_check( self, dummy: Dict ) -> Dict:
        kase = Case(dummy)
        # Check what it needs
        for itag in self.REQUIRED_FIELDS:
            if kase[itag] is None:
                raise NodeDataError(f'Field "{itag}" is required')

        # Include what keywords it adds
        kase.data.setdefault('metadata', {}).setdefault('motif_picker', [])
        return kase.data

    def single_execute( self, data: Dict ) -> Dict:
        kase = Case(data)
        result = {'id': f'motif_{self.identifier}'}

        # Create a working folder
        folders = kase.undirected_path.joinpath(result['id'])
        folders.mkdir(parents=True, exist_ok=True)
        result['data_dir'] = str(folders)

        # Load Structure and create eigens
        try:
            motifs = self.Motif(*reverse_motif(self.log, self.source, self.selection, self.attach,
                                          self.hotspot,  self.identifier, self.binder))
        except StructuralError as se:
            raise NodeDataError(str(se))
        result['motifs'] = motifs
        # Attach data and return
        kase.data.setdefault('metadata', {}).setdefault('motif_picker', []).append(result)
        return kase.data
