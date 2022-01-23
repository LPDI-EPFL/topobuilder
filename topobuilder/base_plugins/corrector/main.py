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
import copy
# import sys

# External Libraries

# This Library
from topobuilder.case import Case
from topobuilder.workflow import Node, NodeDataError


__all__ = ['corrector']


class corrector( Node ):
    """Applies corrections to the placements of the secondary structures in a :term:`FORM`.

    This affects on placement and angles of the different secondary structures. Corrections are defined through
    a controlled vocabulary.

    .. caution::
        Structural motifs imported through the :class:`.Node` :class:`.motif_picker` will impose constraints on
        the way structures can move with respect to each other, thus making some corrections impossible to
        fulfill.

    .. admonition:: To Developers

        Whenever a plugin has to apply geometric and positional changes to the :term:`FORM`, it should be
        done through this class.

    :param corrections: Per secondary structure or per layer corrections to be applied.

    :raises:
        :NodeDataError: On **check**. If the required fields to be executed are not there.
        :NodeDataError: On **execution**. If the requested corrections do not meet the criteria imposed by the
            :term:`FORM`.

    """
    REQUIRED_FIELDS = ('topology.architecture', )
    RETURNED_FIELDS = ()
    VERSION = 'v1.0'

    def __init__( self, tag: int,
                  corrections: Optional[Union[str, Dict, Path, List[Union[str, Path]]]] = None ):
        super(corrector, self).__init__(tag)

        # Make sure we have a list of corrections.
        self.corrections = corrections
        if self.corrections is not None:
            if not isinstance(self.corrections, list):
                self.corrections = [self.corrections, ]
            for i, c in enumerate(self.corrections):
                if isinstance(c, str):
                    self.corrections[i] = Path(c)
        else:
            self.corrections = []

    def single_check( self, dummy: Dict ) -> Dict:
        kase = Case(dummy)

        # Check what it needs
        for itag in self.REQUIRED_FIELDS:
            if kase[itag] is None:
                raise NodeDataError(f'Field "{itag}" is required')

        # Include what keywords it adds (in this instance, nothing)
        return kase.data

    def single_execute( self, data: Dict ) -> Dict:
        kase = Case(data)

        crr = copy.deepcopy(self.corrections)

        # See if there are extra corrections attached to the case itself
        krr = kase['metadata.corrections']
        krr = [] if krr is None else krr
        crr.extend(krr)

        # Apply each set of corrections
        for c in crr:
            self.log.info('Applying correction: {0}\n'.format(c))
            kase = kase.apply_corrections(c)

        self.log.debug(f'Applied a total of {len(crr)} corrections.')
        self.log.debug(f'{len(krr)} from within the Case definition.')
        self.log.debug(f'{len(self.corrections)} from protocol-provided data.')
        return kase.data
