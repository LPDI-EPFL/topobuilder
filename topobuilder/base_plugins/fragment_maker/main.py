# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>
.. codeauthor:: Zander Harteveld <zandermilanh@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from pathlib import Path
from typing import Optional, Tuple, List, Union, Dict
import sys

# External Libraries
import pandas as pd

# This Library
from topobuilder.workflow import Node, NodeDataError, NodeMissingError
from topobuilder.case import Case
import topobuilder.core as TBcore
import topobuilder.utils as TButil
from rstoolbox.io import parse_rosetta_fragments, write_rosetta_fragments, write_fragment_sequence_profiles


__all__ = ['fragment_maker']


class fragment_maker( Node ):
    """Creates or mixes fragments that are needed in multiple Rosetta protocols. Mutliple ways of creating fragments
    are possible through different protocols.

    .. note::
        Currently, solely the ``loop_fragment`` protocol is implemented.

    .. caution::
        In order to create fragments with the ``loop_fragment`` protocol, the :mod:`.loop_fragments` plugin
        needs to be set in the :class:`.Pipeline`.

    :param protocol: Fragment creation protocol to be used.
    :param script: Rosetta script to pick fragments.

    :raises:
        :NodeDataError: On **check**. If the required fields to be executed are not there.
        :NodeDataError: On **execution**. If the :class:`.Case` contains anything other than one defined connectivity.
        :NodeMissingError: On **exection**. If required variable inputs are not there.
    """
    REQUIRED_FIELDS = ('metadata.loop_fragments', 'metadata.loop_lengths')
    RETURNED_FIELDS = ('metadata.fragments')
    VERSION = 'v1.0'

    def __init__( self, tag: int,
                  protocol: str,
                  script: Optional[Union[Path, str]] = None):
        super(fragment_maker, self).__init__(tag)

        self.protocol = protocol
        self.script   = script

    def single_check( self, dummy: Dict ) -> Dict:
        kase = Case(dummy)

        # Check what it needs
        for itag in self.REQUIRED_FIELDS:
            if kase[itag] is None:
                raise NodeDataError(f'Field "{itag}" is required')

        # Include what keywords it adds (in this instance, nothing)
        kase.data.setdefault('metadata', {}).setdefault('fragments', {})
        return kase.data

    def single_execute( self, data: Dict ) -> Dict:
        case = Case(data)
        data = {'protocol': self.protocol, 'files': []}

        # Fragments can only be made for a full, reoriented Case.
        if case.connectivity_count > 1:
            raise NodeDataError('FunFolDes can only be applied to one connectivity.')

        # Generate the folder tree for a single connectivity.
        folders = case.connectivities_paths[0].joinpath('fragment_maker')
        folders.mkdir(parents=True, exist_ok=True)
        checkpoint = folders.joinpath('checkpoint.json')

        # Check if checkpoint exists, retrieve and skip
        reload = TButil.checkpoint_in(self.log, checkpoint)
        if reload is not None and reload['protocol'] == self.protocol:
            case.data['metadata']['fragments'] = reload
            return case

        # Switch depending on the fragment_protocol
        if self.protocol == 'loop_master':
            frags3, frags9, profile = self.loop_master_protocol(case, folders)
            data['files'] = (frags3, frags9)
            data['profile'] = profile
        if self.protocol == 'loopgroup_master':
            #data['files'] = self.loopgroup_master_protocol(case, folders)
            frags3, frags9, profile = self.loop_master_protocol(case, folders)
            data['files'] = (frags3, frags9)
            data['profile'] = profile

        # Store data
        case.data['metadata']['fragments'] = data

        # Checkpoint save
        TButil.checkpoint_out(self.log, checkpoint, data)
        return case


    def loop_master_protocol( self, case: Case, folders: Path ) -> Tuple[str, str]:
        """
        """
        lf = case['metadata.loop_fragments']
        if lf is None:
            raise NodeMissingError('Data that should be loaded through loop_master is not found.')

        for i, loop in enumerate(lf):
            if i == 0:
                ff3 = parse_rosetta_fragments(loop['fragfiles'][0])
                ff9 = parse_rosetta_fragments(loop['fragfiles'][1])
                df3 = [pd.read_csv(str(loop['fragfiles'][0]) + '.csv'), ]
                df9 = [pd.read_csv(str(loop['fragfiles'][1]) + '.csv'), ]
            else:
                df3.append(pd.read_csv(str(loop['fragfiles'][0]) + '.csv'))
                df9.append(pd.read_csv(str(loop['fragfiles'][1]) + '.csv'))
                ff3 = ff3.add_fragments(parse_rosetta_fragments(loop['fragfiles'][0]), ini=int(loop['edges']['ini']), how='append')
                ff9 = ff9.add_fragments(parse_rosetta_fragments(loop['fragfiles'][1]), ini=int(loop['edges']['ini']), how='append')

        TButil.plot_fragment_templates(self.log, pd.concat(df3), pd.concat(df9), folders.joinpath('template_fragment_profile'))

        small_file = write_rosetta_fragments(ff3.top_limit(lf[-1]['edges']['end']), prefix=folders.joinpath('small'), strict=True)
        self.log.info(f'Writing small fragment file: {small_file}')
        large_file = write_rosetta_fragments(ff9.top_limit(lf[-1]['edges']['end']), prefix=folders.joinpath('large'), strict=True)
        self.log.info(f'Writing large fragment files: {large_file}')

        # Create fragment design profile
        profile_file = folders.joinpath('frag_profile.pssm')
        dfs9all = parse_rosetta_fragments(str(folders)+"/large.200.9mers")
        write_fragment_sequence_profiles(dfs9all, filename=profile_file)
        self.log.info(f'Writing sequence profile: {str(folders)+"/large.200.9mers"}')

        return small_file, large_file, profile_file

    def loopgroup_master_protocol( self, case: Case, folders: Path ) -> Tuple[str, str]:
        """
        """
        lf = case['metadata.loop_fragments']
        if lf is None:
            raise NodeMissingError('Data that should be loaded through loop_master is not found.')
        for i, loop in enumerate(lf):
            if i == 0:
                ff3 = parse_rosetta_fragments(loop['fragfiles'][0])
                ff9 = parse_rosetta_fragments(loop['fragfiles'][1])
                df3 = [pd.read_csv(str(loop['fragfiles'][0]) + '.csv'), ]
                df9 = [pd.read_csv(str(loop['fragfiles'][1]) + '.csv'), ]
            else:
                df3.append(pd.read_csv(str(loop['fragfiles'][0]) + '.csv'))
                df9.append(pd.read_csv(str(loop['fragfiles'][1]) + '.csv'))
                ff3 = ff3.add_fragments(parse_rosetta_fragments(loop['fragfiles'][0]), ini=int(loop['edges']['ini']), how='append')
                ff9 = ff9.add_fragments(parse_rosetta_fragments(loop['fragfiles'][1]), ini=int(loop['edges']['ini']), how='append')

        TButil.plot_fragment_templates(self.log, pd.concat(df3), pd.concat(df9), folders.joinpath('template_fragment_profile'))

        small_file = write_rosetta_fragments(ff3.top_limit(lf[-1]['edges']['end'] - 3), prefix=folders.joinpath('small'), strict=True)
        self.log.info(f'Writing small fragment file: {small_file}\n')
        large_file = write_rosetta_fragments(ff9.top_limit(lf[-1]['edges']['end'] - 3), prefix=folders.joinpath('large'), strict=True)
        self.log.info(f'Writing large fragment files: {large_file}\n')

        # Create fragment design profile
        profile_file = folders.joinpath('frag_profile.pssm')
        dfs9all = parse_rosetta_fragments(str(folders)+"/large.200.9mers")
        write_fragment_sequence_profiles(dfs9all, filename=profile_file)
        self.log.info(f'Writing sequence profile: {str(folders)+"/large.200.9mers"}')

        return small_file, large_file, profile_file
