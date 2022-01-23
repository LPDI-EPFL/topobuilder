# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>
.. codeauthor:: Zander Harteveld <zandermilanh@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import List, Dict, Optional
import sys

# External Libraries
import numpy as np
import pandas as pd
try:
    from SBI.structure import PDB, Frame3D
    import SBI.core as SBIcr
except ImportError:
    class Frame3D():
        pass

# This Library
from topobuilder.workflow import Node, NodeDataError
from topobuilder.case import Case
import topobuilder.core as TBcore
import topobuilder.utils as TButil
from .parametric import SSEArchitect


__all__ = ['builder']


class builder( Node ):
    """Builds and adds a :term:`SKETCH` from given :term:`FORM` string to the :class:`.Case` using ideal SSE elements.

    If corrections are available and specified, these will be applied onto the :term:`SKETCH`.

    .. caution::
        In order to apply secondary structure or per layer corrections, the :mod:`.corrector` plugin
        needs to be set in the :class:`.Pipeline`.

    :param connectivity: Expected secondary structure connectivity. *Important*: at the moment only a single
                         connectivity supported (default: True).
    :param motif: Expected Motif to be added to the :term:`SKETCH` (default: False).
    :param pick_aa: Desired amino acid type to use for the :term:`SKETCH` sequence. If not specified, it will
                    use pseudorandomly assign amino acid types based on secondary structure propensity scores.
    :param write2disc: Dump the :term:`SKETCH` (default: True).

    :raises:
        :NodeDataError: On **check**. If the required fields to be executed are not there.
        :NodeDataError: On **execution**. If the :class:`.Case` contains anything other than one defined connectivity.
    """
    REQUIRED_FIELDS = ('topology.architecture', 'topology.connectivity')
    RETURNED_FIELDS = ()
    VERSION = 'v1.0'

    def __init__( self, tag: int,
                  connectivity: Optional[bool] = True,
                  motif: Optional[bool] = False,
                  pick_aa: Optional[str] = None,
                  write2disc: Optional[str] = True):
        super(builder, self).__init__(tag)

        self.connectivity = connectivity
        self.motif = motif
        self.pick_aa = pick_aa
        self.write2disc = write2disc

    def single_check( self, dummy: Dict ) -> Dict:
        kase = Case(dummy)

        # Check what it needs
        for itag in self.REQUIRED_FIELDS:
            if kase[itag] is None:
                raise NodeDataError(f'Field "{itag}" is required')

        # Include what keywords it adds (in this instance, nothing)
        # Here Nothing
        return kase.data

    def single_execute( self, data: Dict ) -> Dict:
        case = Case(data)

        # Apply connectivity?
        if self.connectivity:
            if case.connectivity_count == 0:
                raise NodeDataError('Minimum a single connectivity must be provided.')
            if case.connectivity_count > 1:
                raise NodeDataError('Only single connectivity cases can be build.')
            case = case.cast_absolute().apply_topologies()[0]
            ofile = case.connectivities_paths[0].joinpath('directed_sketch.pdb')
            bfile = case.connectivities_paths[0].joinpath('binder.pdb')
            ifile = case.connectivities_paths[0].joinpath('input.pdb')
        else:
            case = case.cast_absolute()
            ofile = case.main_path.joinpath('architecture').joinpath('undirected_sketch.pdb')

        # for i, j, sse in case:
        #     if 'atoms' in sse['metadata'] and not TBcore.get_option('system', 'overwrite'):
        #         self.log.debug(f'{case.name}.{sse["id"]} already has atoms defined\n')
        #         continue
        #
        #     self.log.debug(f'Building coordinates for {case.name}.{sse["id"]}\n')
        #     case.data['topology']['architecture'][i][j] = self.make_structure(sse, self.pick_aa)

        # Insert motif?
        if self.motif:
            if not case.data['metadata']['motif_picker']:
                raise NodeDataError('Motif must be provided and pre-picked through the motif_picker.')
            else:
                # Include what keywords it adds (in this instance, nothing)
                case.data.setdefault('metadata', {}).setdefault('binder', {})

                for n, motif_data in enumerate(case.data['metadata']['motif_picker']):
                    motif, binder, hotspots, attach, selection, identifier = motif_data['motifs']
                    motif = pd.DataFrame(motif, columns=['auth_comp_id', 'auth_atom_id', 'auth_seq_id', 'auth_asym_id',
                                                         'sse_id', 'internal_num', 'Cartn_x', 'Cartn_y', 'Cartn_z']
                                                         ).reset_index(drop=True)
                    binder = pd.DataFrame(binder, columns=['auth_comp_id', 'auth_atom_id', 'auth_seq_id', 'auth_asym_id',
                                                           'Cartn_x', 'Cartn_y', 'Cartn_z']).reset_index(drop=True)
                    self.log.debug(f'processing motif id: {identifier}')

                    mcolumns, bcolumns = motif.columns.tolist(), binder.columns.tolist()
                    initial = PDB( pd.concat([motif[bcolumns], binder[bcolumns]], sort=False) )

                    # Get segments of interest
                    segments = []
                    for i,j,sse in case:
                        if sse['id'] in attach:
                            segment = pd.DataFrame( sse['metadata']['atoms'],
                                                    columns=['residue', 'auth_atom_id', 'resi_id',
                                                             'Cartn_x', 'Cartn_y', 'Cartn_z'] )
                            segment = segment.assign(sse_id=[sse['id']]*len(segment))
                            segments.append(segment)
                    segments = pd.concat(segments, sort=False)

                    # Prepare for alignment
                    att_segments = segments.drop(['Cartn_x', 'Cartn_y', 'Cartn_z'], axis=1)
                    att_motifs   = motif.drop(['Cartn_x', 'Cartn_y', 'Cartn_z'], axis=1)
                    fl_segments  = segments[['Cartn_x', 'Cartn_y', 'Cartn_z']].values
                    fl_motif     = motif[['Cartn_x', 'Cartn_y', 'Cartn_z']].values
                    ca_segments  = segments[segments['auth_atom_id'] == 'CA'][['Cartn_x', 'Cartn_y', 'Cartn_z']].values
                    ca_motifs    = motif[motif['auth_atom_id'] == 'CA'][['Cartn_x', 'Cartn_y', 'Cartn_z']].values

                    # Align
                    self.log.debug(f'Aligning motif of length {len(ca_segments)} onto segement of length {len(ca_motifs)}')
                    fl_motif_ali = self.align_by_Kabsch(fl_motif, ca_segments, moving_guide=ca_motifs)

                    # Recreate frame
                    fl_motif_ali = pd.concat([pd.DataFrame(att_motifs, columns=att_motifs.columns.tolist()).reset_index(drop=True),
                                              pd.DataFrame(fl_motif_ali, columns=['Cartn_x', 'Cartn_y', 'Cartn_z']).reset_index(drop=True)],
                                              sort=False, axis=1)

                    # Correct in case
                    for i, j, sse in case:
                        if 'atoms' in sse['metadata'] and not TBcore.get_option('system', 'overwrite'):
                            self.log.debug(f'{case.name}.{sse["id"]} already has atoms defined')
                            continue
                        if sse['id'] in attach: # part of the motif
                            self.log.debug(f'Writing motif coordinates for {case.name}.{sse["id"]}')
                            fl_motif_ali_part = fl_motif_ali[fl_motif_ali['sse_id'] == sse['id']]
                            motif_list = list(fl_motif_ali_part[['auth_comp_id', 'auth_atom_id', 'auth_seq_id', #'internal_num',
                                                                 'Cartn_x', 'Cartn_y', 'Cartn_z']].values)
                            case.data['topology']['architecture'][i][j] = self.make_motif(sse, motif_list)
                        else: # not motif
                            self.log.debug(f'Writing coordinates for {case.name}.{sse["id"]}')
                            case.data['topology']['architecture'][i][j] = self.make_structure(sse, self.pick_aa )

                    if binder is not None:
                        self.log.debug(f'Aligning binder next to motif')
                        att_binder    = binder.drop(['Cartn_x', 'Cartn_y', 'Cartn_z'], axis=1)
                        fl_binder     = binder[['Cartn_x', 'Cartn_y', 'Cartn_z']].values
                        fl_binder_ali = self.align_by_Kabsch(fl_binder, ca_segments, moving_guide=ca_motifs)
                        fl_binder_ali = pd.concat([pd.DataFrame(att_binder, columns=att_binder.columns.tolist()).reset_index(drop=True),
                                                   pd.DataFrame(fl_binder_ali, columns=['Cartn_x', 'Cartn_y', 'Cartn_z']).reset_index(drop=True)],
                                                   sort=False, axis=1)
                        # Chain id guided by *auth_asym_id*
                        fl_binder_ali = PDB( fl_binder_ali[['auth_comp_id', 'auth_atom_id', 'auth_asym_id', 'auth_seq_id',
                                                            'Cartn_x', 'Cartn_y', 'Cartn_z']] )
                        #for i, j, sse in case:
                        #    case.data['topology']['architecture'][i][j] = self.add_binder(sse, list(fl_binder_ali[bcolumns].values) )
                        #    break

                        case.data['metadata']['binder'][str(identifier)] = fl_binder_ali
        else:
            for i, j, sse in case:
                if 'atoms' in sse['metadata'] and not TBcore.get_option('system', 'overwrite'):
                    self.log.debug(f'{case.name}.{sse["id"]} already has atoms defined')
                    continue

                self.log.debug(f'Building coordinates for {case.name}.{sse["id"]}')
                case.data['topology']['architecture'][i][j] = self.make_structure(sse, self.pick_aa)

        if self.write2disc:
            ofile.parent.mkdir(parents=True, exist_ok=True)
            structure, _ = TButil.build_pdb_object( self.log, case.ordered_structures, 2 )

            self.log.notice(f'Writing structure {ofile}')
            structure.write(output_file=str(ofile), format='pdb', clean=True,
                            force=TBcore.get_option('system', 'overwrite'))

            if self.motif and binder is not None:
                bfile.parent.mkdir(parents=True, exist_ok=True)
                self.log.notice(f'Writing binder {bfile}')
                fl_binder_ali.write(output_file=str(bfile), format='pdb', clean=True,
                                    force=TBcore.get_option('system', 'overwrite'))

                ifile.parent.mkdir(parents=True, exist_ok=True)
                self.log.notice(f'Writing input {ifile}')
                initial.write(output_file=str(ifile), format='pdb', clean=True,
                                    force=TBcore.get_option('system', 'overwrite'))

        return case


    def make_structure( self, sse: Dict, pick_aa: Optional[str] = None ) -> Case:
        """
        """
        structure = SSEArchitect(sse, type=sse['type'], pick_aa=pick_aa).pdb
        sse['metadata'].setdefault('atoms', None)
        sse['metadata']['atoms'] = list(structure[['auth_comp_id', 'auth_atom_id', 'auth_seq_id',
                                                   'Cartn_x', 'Cartn_y', 'Cartn_z']].values)
        return sse


    def make_motif( self, sse: Dict, motif: List[str]) -> Case:
        """
        """
        sse['metadata'].setdefault('atoms', None)
        sse['metadata']['atoms'] = motif
        return sse


    def add_binder( self, sse: Dict, binder: List[str]) -> Case:
        """
        """
        sse['metadata'].setdefault('binder', None)
        sse['metadata']['binder'] = binder
        return sse


    def align_by_Kabsch(self, moving_selection, target_guide, moving_guide = None):
        """
        """
        if moving_guide is None: moving_guide = np.copy(moving_selection)
        assert len(target_guide) == len(moving_guide)

        guide_center  = np.mean(moving_guide, axis = 0, dtype = np.float64)
        moving_center = np.mean(moving_selection, axis = 0, dtype = np.float64)
        target_center = np.mean(target_guide, axis = 0, dtype = np.float64)

        M = self.kabsch(moving_guide - guide_center, np.copy(target_guide) - target_center)

        moving_selection -= guide_center
        moving_selection = np.dot(moving_selection, M)
        moving_selection += target_center

        return moving_selection


    def kabsch(self, P, Q):
        """
        The optimal rotation matrix U is calculated and then used to rotate matrix
        P unto matrix Q so the minimum root-mean-square deviation (RMSD) can be
        calculated.
        Using the Kabsch algorithm with two sets of paired point P and Q,
        centered around the center-of-mass.
        Each vector set is represented as an NxD matrix, where D is the
        the dimension of the space.
        The algorithm works in three steps:
        - a translation of P and Q
        - the computation of a covariance matrix C
        - computation of the optimal rotation matrix U
        http://en.wikipedia.org/wiki/Kabsch_algorithm
        Parameters:
        P -- (N, number of points)x(D, dimension) matrix
        Q -- (N, number of points)x(D, dimension) matrix
        Returns:
        U -- Rotation matrix
        """

        # Computation of the covariance matrix
        C = np.dot(np.transpose(P), Q)

        # Computation of the optimal rotation matrix
        # This can be done using singular value decomposition (SVD)
        # Getting the sign of the det(V)*(W) to decide
        # whether we need to correct our rotation matrix to ensure a
        # right-handed coordinate system.
        # And finally calculating the optimal rotation matrix U
        # see http://en.wikipedia.org/wiki/Kabsch_algorithm
        V, S, W = np.linalg.svd(C)
        d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

        if d:
            S[-1] = -S[-1]
            V[:, -1] = -V[:, -1]

        # Create Rotation matrix U
        U = np.dot(V, W)
        return U
