# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>
.. codeauthor:: Zander Harteveld <zandermilanh@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import os
from pathlib import Path
from typing import Optional, List, Dict
import sys
# import textwrap

# External Libraries

# This Library
from topobuilder.workflow import Node, NodeDataError, NodeMissingError
from topobuilder.case import Case
import topobuilder.core as TBcore
import topobuilder.utils as TButil
from . import utils


__all__ = ['hybridize']


class hybridize( Node ):
    """Run Rosettas hybridize protocol to generate designs.

    .. caution::
        Due to the ``FastDesignMover``, this :class:`funfoldes` may take *a lot* of time. If possible, please
        use the ``slurm.use`` configuration. In case this is not possible, you may reduce the number of decoys
        to be generated via the `nstruct` parameter option.


    :param nstruct: Number of decoys to be generated (default: 2000).
    :param natbias: Score function bias towards per secondary structure types (default: 2.5).
    :param layer_design: If :class:`funfoldes` should a layer design approach (default: True).

    :raises:
        :NodeDataError: On **check**. If the required fields to be executed are not there.
        :NodeMissingError: On **exection**. If required variable inputs are not there.
    """
    REQUIRED_FIELDS = ('metadata.fragments','metadata.loop_lengths')
    RETURNED_FIELDS = ('metadata.hybridize')
    VERSION = 'v1.0'

    def __init__( self, tag: int,
                  nstruct: Optional[int] = 2000,
                  natbias: Optional[float] = 2.5,
                  layer_design: Optional[bool] = True) -> Case:
        super(hybridize, self).__init__(tag)

        self.nstruct = nstruct
        self.natbias = natbias
        self.layer_design = layer_design

    def single_check( self, dummy: Dict ) -> Dict:
        kase = Case(dummy)

        # Check what it needs
        for itag in self.REQUIRED_FIELDS:
            if kase[itag] is None:
                raise NodeDataError(f'Field "{itag}" is required')

        # Include what keywords it adds (in this instance, nothing)
        kase.data.setdefault('metadata', {}).setdefault('hybridize', {})
        return kase.data

    def single_execute( self, data: Dict ) -> Dict:
        case = Case(data)

        data = {'script': {'assembly': '', 'design': ''},
                'cmd':
                    {'assembly': [Path(TBcore.get_option('rosetta', 'scripts')), '-parser:protocol'],
                     'design': [Path(TBcore.get_option('rosetta', 'scripts')), '-parser:protocol']},
                'silent_files': {'assembly': [], 'design': []},
                'minisilent': {'assembly': '', 'design': ''}
                }

        # Generate the folder tree for a single connectivity.
        wpaths = utils.folder_structure(case)

        # Check if checkpoint exists, retrieve and skip
        reload = TButil.checkpoint_in(self.log, wpaths['checkpoint'])
        if reload is not None:
            case.data['metadata']['hybridize'] = reload
            return case

        # We need to check that the rosetta_scripts executable is available
        if not data['cmd']['folding'][0].is_file() or not os.access(str(data['cmd']['assembly'][0]), os.X_OK):
            raise NodeMissingError(f'Cannot find executable {data["cmd"]["assembly"][0]}')

        # Build the structure
        utils.build_template_sketch( self.log, case, wpaths['pdb'] )

        # Make the Folding and Design RScripts
        data = utils.make_scripts(self.log, case, wpaths, data, self.natbias, self.layer_design)

        # Finish command
        data = utils.commands(case, self.nstruct, data, wpaths)

        # Execute
        data = utils.execute(self.log, data, wpaths)

        # Update metadata
        data = utils.update_data(self.log, data, wpaths)

        # Checkpoint save
        TButil.checkpoint_out(self.log, wpaths['checkpoint'], data)
        case.data['metadata']['funfoldes'] = data

        return case
