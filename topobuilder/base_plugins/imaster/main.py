# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>
.. codeauthor:: Zander Harteveld <zandermilanh@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import Optional, List, Dict, Set
from string import ascii_uppercase
from operator import itemgetter
from subprocess import run, DEVNULL
from pathlib import Path
from itertools import cycle
import shlex
import math
import os
import sys
import glob

# External Libraries
import numpy as np
import pandas as pd
import scipy as sc
import networkx as nx

# This Library
from topobuilder.case import Case
from topobuilder.workflow import Node, NodeOptionsError, NodeDataError
import topobuilder.core as TBcore
import topobuilder.utils as TButil
from topobuilder.workflow import TBplugins
from .analysis import get_steps


__all__ = ['imaster']


class imaster( Node ):
    """Searches for defined per secondary structure or layer corrections by matching layer-wise against
    a defined :ref:`MASTER` database. If a beta-layer is present, this will be the starting layer to be
    corrected. Multiple beta-layers are sorted by size (number of strands) and the largest one is set as
    starting layer. From there on, all next layers will be corrected paire-wise, e.g. the previous
    corrected layer together with the next layer.

    .. note::
        Multiple correction shemes can be set. The default correction sheme is layer tilt around the x-axis and
        layer points (shifts) along the y-axis.

    .. caution::
        Depending on the size of the queried :ref:`MASTER` database, this may take *a lot of time*.
        If possible, please use the ``slurm.use`` configuration. In case this is not possible, you may
        reduce the size of the database by setting ``master.pds`` pointing to fewer structures to search over.

    .. admonition:: To Developers

        Due to the possibilty of external :class:`.Node`, main function is located in the :mod:`.imaster` module.

    :param rmsd: RMSD threshold for master search (default: 5.0).
    :param bin: Starting bin for master corrections (default: mid).
    :param step: Path of how the layers are searched for corrections.
    :param subsampling: Subsampling matches to accelerated calculations.
    :param corrections: File containing the corrections.
    :param mirror_beta_twist: Invert beta sheet twist (default: False).
    :param mirror_beta_shear: Invert beta sheet shear between layers (default: False).
    :param correction_type: Analysis method ['network', 'mode'] to use for corrections (default: 'network').
    :param correction_check: Check if corrections have been correctly applied across the layers (default: True).

    :raises:
        :NodeOptionsError: On **initialization**. If a reserved key is provided as a subname.
        :NodeDataError: On **check**. If the required fields to be executed are not there.

    """
    RESERVED_KEYWORDS = ['imaster']
    REQUIRED_FIELDS = ('configuration.name', 'topology.connectivity')
    RETURNED_FIELDS = ('metadata.corrections', 'metadata.imaster', 'metadata.bin')
    VERSION = 'v1.0'

    def __init__( self, tag: int,
                  correctives:  Optional[Dict] = {'rotate': ['x',], 'translate': ['y',]},
                  rmsd: Optional[float] = 5.0,
                  bin: Optional[str] = 'mid',
                  step: Optional[int] = None,
                  subsampling: Optional[int] = None,
                  mirror_beta_twist: Optional[bool] = False,
                  mirror_beta_shear: Optional[bool] = False,
                  corrections: Optional[Dict] = dict(),
                  correction_type: Optional[str] = 'network',
                  correction_check: Optional[bool] = True ):
        super(imaster, self).__init__(tag)

        #self.cases = cases
        self.correctives = correctives
        self.rmsd = rmsd
        self.bin = bin
        self.step = step
        self.subsampling = subsampling
        self.mirror_beta_twist = mirror_beta_twist
        self.mirror_beta_shear = mirror_beta_shear
        self.corrections = corrections if TBcore.get_option('system', 'jupyter') else {}
        self.correction_type = correction_type
        self.correction_check = correction_check


    def single_check( self, dummy: Dict ) -> Dict:
        kase = Case(dummy)

        # Check what it needs
        for itag in self.REQUIRED_FIELDS:
            if kase[itag] is None:
                raise NodeDataError(f'Field "{itag}" is required')

        # Include what keywords it adds (in this instance, nothing)
        kase.data.setdefault('metadata', {}).setdefault('imaster', {})
        kase.data.setdefault('metadata', {}).setdefault('corrections', [])
        kase.data.setdefault('metadata', {}).setdefault('bin', self.bin)
        return kase.data


    def single_execute( self, data: Dict ) -> Dict:
        kase = Case(data)

        # Include what keywords it adds (in this instance, nothing)
        kase.data.setdefault('metadata', {}).setdefault('imaster', {})
        kase.data.setdefault('metadata', {}).setdefault('corrections', [])
        kase.data.setdefault('metadata', {}).setdefault('bin', self.bin)

        # Interactive MASTER smoothing is only applied to a Case with one single connectivity
        if kase.connectivity_count > 1:
            raise NodeDataError('Interactive MASTER smoothing can only be applied to one connectivity.')
        if not kase['configuration.reoriented'] and kase.connectivity_count == 1:
            kase = kase.cast_absolute().apply_topologies()[0]

        # Generate the folder tree for a single connectivity.
        wfolder = kase.connectivities_paths[0].joinpath('imaster')
        wfolder.mkdir(parents=True, exist_ok=True)
        current_case_file = kase.cast_absolute().write(wfolder.joinpath('current'))

        # Find steps: Tops we will submit 2-layer searches
        steps = get_steps([x[-1] == 'E' for x in kase.architecture_str.split('.')])
        steps = [steps[self.step], ] if self.step is not None and TBcore.get_option('system', 'jupyter') else steps

        # Work by layers
        data_prev, corr_prev = None, None
        done_l = set()
        for i, step in enumerate(steps):
            # Step working directory
            stepfolder = wfolder.joinpath('step{:02d}'.format(i + 1))
            stepfolder.mkdir(parents=True, exist_ok=True)
            query = stepfolder.joinpath('imaster.query{:02d}.pdb'.format(i + 1))
            checkpoint = stepfolder.joinpath('checkpoint.json')

            reload = TButil.checkpoint_in(self.log, checkpoint)
            if reload is not None:
                kase.data['metadata']['imaster'].setdefault('step{:02d}'.format(i + 1), reload)
                kase.data['metadata']['corrections'].append(reload['corrections'])
                self.corrections.update(reload['corrections'])
                done_l.update(reload['layers'])
                data_prev = pd.read_csv(stepfolder.joinpath('geometry.csv'))
                # CKase = CKase.apply_corrections(corrections)
                continue

            # Apply corrections from previous steps and rebuild
            CKase = Case(kase).apply_corrections(self.corrections)
            with TBcore.on_option_value('system', 'overwrite', True):
                node = getattr(TBplugins.source.load_plugin('builder'), 'builder', None)(connectivity=True, tag=0)
                CKase = Case(node.single_execute(CKase.data))

            # Generate structure query and get layer displacements
            layers = set(itemgetter(*step)(ascii_uppercase))
            sses = [sse for sse in CKase.ordered_structures if sse['id'][0] in layers]
            structure, cends = TButil.build_pdb_object( self.log, sses, 3)
            self.log.notice(f'New file: Writing structure {query}')
            structure.write(output_file=str(query), format='pdb', clean=True, force=True)

            flip = cycle([CKase['configuration.flip_first'], not CKase['configuration.flip_first']])
            counts = np.asarray([sse['length'] for sse in CKase.ordered_structures])
            cends = np.cumsum(counts)
            cstrs = cends - counts + 1

            rules = list(zip([sse['id'] for sse in CKase.ordered_structures],
                             list(zip(cstrs, cends)),
                             list(next(flip) for _ in range(len(CKase.ordered_structures)))))
            extras = TButil.pdb_geometry_from_rules(query, rules, self.log)

            # MASTER search
            createpds = TButil.createPDS(query)
            self.log.notice(f'EXECUTE: {" ".join([str(x) for x in createpds])}')
            run(createpds, stdout=DEVNULL)
            masters = TButil.master_best_each(self.log, query.with_suffix('.pds'), stepfolder.joinpath('_master'), self.rmsd)
            data = self.submit_searches(masters, stepfolder, current_case_file, '.'.join([x['id'] for x in sses]))
            data = self.calc_corrections(data, kase, set(data['layers']), done_l, extras, rules, bin=self.bin, data_prev=data_prev, corr_prev=corr_prev)
            corr_prev = data['corrections']
            data_prev = pd.read_csv(stepfolder.joinpath('geometry.csv'))

            kase.data['metadata']['imaster'].setdefault('step{:02d}'.format(i + 1), data)
            TButil.checkpoint_out(self.log, checkpoint, data)
            kase.data['metadata']['corrections'].append(data['corrections'])
            done_l.update(data['layers'])
            self.corrections.update(data['corrections'])

        # Correct beta twist
        if self.mirror_beta_twist is True:
            self.log.debug('Mirroring beta tilt corrections')
            for dd in kase.data['metadata']['corrections']:
                for sse in dd.keys():
                    dd[sse]['tilt']['x'] *= -1.

        if self.mirror_beta_shear is True:
            self.log.debug('Mirroring beta shear corrections')
            for dd in kase.data['metadata']['corrections']:
                for sse in dd.keys():
                    try:
                        dd[sse]['tilt']['z'] *= -1.
                    except:
                        continue

        return kase


    def submit_searches( self, cmd: List[str], wdir: Path, current_case_file: Path, current_sse: str ) -> Dict:
        """
        """
        unimaster = wdir.joinpath('match.master')
        imaster = Path(__file__).parent.joinpath('imaster.py')
        unidata = wdir.joinpath('geometry.csv')
        if unimaster.is_file() and unidata.is_file():
            return {'matches': unimaster, 'stats': unidata, 'corrections': None,
                    'layers': list(set([x[0] for x in current_sse.split('.')]))}
        if not TBcore.get_option('slurm', 'use'):
            self.no_slurm(cmd, current_case_file, current_sse, unimaster, imaster,
                          unidata.with_suffix(''))
        else:
            self.with_slurm(cmd, current_case_file, current_sse, unimaster, imaster, unidata)

        return {'matches': unimaster, 'stats': unidata, 'corrections': None,
                'layers': list(set([x[0] for x in current_sse.split('.')]))}


    def calc_corrections( self, data: Dict, case: Case, qlayers: Set, dlayers: Set,
                          extras: pd.DataFrame, rules: list, bin: Optional[str] = 'mid',
                          data_prev: Optional[Dict] = None, corr_prev: Optional[Dict] = None) -> Dict:
        """
        """
        tocorrect = qlayers.difference(dlayers)
        toreference = qlayers.difference(tocorrect)
        if len(tocorrect) > 1:
            raise ValueError('Layers are corrected one by one.')
        tocorrect = list(tocorrect)[0]

        if len(toreference) == 1:
            toreference = list(toreference)[0]
        elif len(toreference) == 0:
            toreference = None
        else:
            for x in toreference:
                found = False
                if abs(ascii_uppercase.find(x) - ascii_uppercase.find(tocorrect)) == 1:
                    toreference = x
                    found = True
                    break
                if not found:
                    toreference = None

        # Load Data, bin, show and addapt if no matches for the given bin.
        bins = ["close", "mid", "far", "extreme"]
        df = pd.read_csv(data['stats'])

        # Correct based on the previous orientations.
        if self.correction_check == True:
            if toreference is not None:
                df_prev = data_prev[data_prev.layer == toreference]
                dfss    = df[(df.sse.isin(df_prev.sse)) & (df.layer.isin(df_prev.layer))]
                for a in ['angles_layer', 'angles_side']:
                    sign_flip = False
                    m_prev = np.mean([abs(r)/r for r in df_prev[a].values])
                    if m_prev != 0.: m_prev = abs(m_prev)/m_prev
                    m_ss   = np.mean([abs(r)/r for r in dfss[a].values])
                    if m_ss != 0.: m_ss = abs(m_ss)/m_ss

                    #if m_ss != 0. and m_prev != 0.:
                    if m_ss != m_prev: sign_flip = True

                    if sign_flip:
                        df[df.layer == tocorrect].loc[:, a] *= -1.
                        self.log.debug(f'{a}: layers current {tocorrect}, previous {toreference} - sign flip: {a}, with value {m_ss} vs. {m_prev}')

        df = df.assign(bin=pd.cut(df['rmsd'], bins=[0, 2, 2.5, 3, 5], labels=bins))
        bin_sizes = df.groupby('bin').count().match.values
        _, _, isBin = TButil.plot_match_bin(self.log, df, Path(data['stats']).parent.joinpath('match_count'),
                                            len(TButil.pds_database(self.log)[1]), ['pdb', 'chain'])

        for i in range(len(bin_sizes)):
            binInt = bins.index(bin)
            if binInt >= len(bins):
                self.log.info('No matches found to use for correction.')
                data['corrections'] = {}
                return data
            if not TBcore.get_option('slurm', 'use'):
                if bin_sizes[binInt] < 1: # avoid taking metrics from very few examples
                    self.log.warning(f'Not enough matches at {bin} ({bin_sizes[binInt]} < 1), trying higher...')
                    binInt += 1
                    bin = bins[binInt]
                    data['bin'] = bin
            elif TBcore.get_option('slurm', 'use'):
                if bin_sizes[binInt] < 1:
                    self.log.warning(f'Not enough matches at {bin} ({bin_sizes[binInt]} < 1), trying higher...')
                    binInt += 1
                    bin = bins[binInt]
                    data['bin'] = bin
            else:
                bin = bins[binInt]
                data['bin'] = bin
                break

        if toreference is None:
            if case.get_type_for_layer(tocorrect) == 'E':
                data['corrections'] = self.first_layer_correction(df, bin, Path(data['stats']).parent)
            else:
                # for now the same as for the betas
                data['corrections'] = self.first_layer_correction(df, bin, Path(data['stats']).parent)
        elif case.get_type_for_layer(toreference) == 'E':
            if case.get_type_for_layer(tocorrect) == 'H':
                if self.correction_type == 'network':
                    self.log.debug('Network correction approach.\n')
                    data['corrections'] = self.nth_layer_correction(df, case, bin, Path(data['stats']).parent, tocorrect, toreference, rules)
                elif self.correction_type == 'mode':
                    self.log.debug('Mode correction approach.\n')
                    data['corrections'], data['prefixes'] = self.alpha_on_beta_correction(df, bin, Path(data['stats']).parent, tocorrect, toreference, case, extras)
                else:
                    raise NodeOptionsError(f'Unknown correction type {self.correction_type}')
            else:
               #data['corrections'], data['prefixes'] = self.beta_on_beta_correction(df, bin, Path(data['stats']).parent, tocorrect, toreference, case, extras, rules)
               if self.correction_type == 'network':
                   self.log.debug('Network correction approach.\n')
                   data['corrections'] = self.nth_layer_correction(df, case, bin, Path(data['stats']).parent, tocorrect, toreference, rules)
               elif self.correction_type == 'mode':
                   self.log.debug('Mode correction approach.\n')
                   data['corrections'], data['prefixes'] = self.beta_on_beta_correction(df, bin, Path(data['stats']).parent, tocorrect, toreference, case, extras, rules)
               else:
                   raise NodeOptionsError(f'Unknown correction type {self.correction_type}')
        elif case.get_type_for_layer(toreference) == 'H':
            if case.get_type_for_layer(tocorrect) == 'H':
                #data['corrections'], data['prefixes'] = self.alpha_on_beta_correction(df, bin, Path(data['stats']).parent, tocorrect, toreference, case, extras)
                if self.correction_type == 'network':
                    self.log.debug('Network correction approach.\n')
                    data['corrections'] = self.nth_layer_correction(df, case, bin, Path(data['stats']).parent, tocorrect, toreference, rules)
                elif self.correction_type == 'mode':
                    self.log.debug('Mode correction approach.\n')
                    data['corrections'], data['prefixes'] = self.alpha_on_beta_correction(df, bin, Path(data['stats']).parent, tocorrect, toreference, case, extras)
                else:
                    raise NodeOptionsError(f'Unknown correction type {self.correction_type}')
            else:
                if self.correction_type == 'network':
                    self.log.debug('Network correction approach.\n')
                    data['corrections'] = self.nth_layer_correction(df, case, bin, Path(data['stats']).parent, tocorrect, toreference, rules)
                elif self.correction_type == 'mode':
                    self.log.debug('Mode correction approach.\n')
                    data['corrections'], data['prefixes'] = self.alpha_on_beta_correction(df, bin, Path(data['stats']).parent, tocorrect, toreference, case, extras)
                else:
                    raise NodeOptionsError(f'Unknown correction type {self.correction_type}')

        if toreference is not None:
            if corr_prev is not None:
                # The case where we have a ab-sandwich, and the helices to adjust
                # we assume that a helix should be rooted by 2 strands.
                if case.get_type_for_layer(toreference) == 'E' and case.get_type_for_layer(tocorrect) == 'H':
                    self.log.info('The consecutive layers are of type EH. Re-checking for correct layer angles automatically.')
                    sses1 = [k for k,v in corr_prev.items()]
                    sses2 = data['corrections'].keys()
                    sse1_cont = []
                    for ssi in range(0, len(sses1) - 1, 2):
                        ssei0, ssei1 = sses1[ssi], sses1[ssi + 1]
                        ssei0_tx  = corr_prev[ssei0]['tilt']['x']
                        ssei1_tx  = corr_prev[ssei1]['tilt']['x']
                        sseim_tx  = (ssei0_tx + ssei1_tx) / 2
                        sse1_cont.append(abs(sseim_tx) / sseim_tx)
                    for ssi, sse2 in enumerate(sses2):
                        self.log.debug(f'{sse2} vs beta')
                        sse2_tx   = data['corrections'][sse2]['tilt']['x']
                        sse2_norm = abs(sse2_tx)/sse2_tx
                        if sse1_cont[ssi] != sse2_norm:
                            self.log.debug(f'Re-correcting {sse2} tilt.')
                            data['corrections'][sse2]['tilt']['x'] *= -1.
                # The beta sandwich case
                if case.get_type_for_layer(toreference) == 'E' and case.get_type_for_layer(tocorrect) == 'E':
                    self.log.info('The 2 consecutive layers are of type E. Re-checking for correct layer angles automatically.')
                    sses1 = corr_prev.keys()
                    sses2 = data['corrections'].keys()
                    for sse1 in sses1:
                        for sse2 in sses2:
                            self.log.debug(f'{sse1} vs {sse2}')
                            if sse1[1] == sse2[1] and sse1[0] != sse2[0]:
                                self.log.debug(f'Comparing {sse1} against {sse2}')
                                sse1_tx   = corr_prev[sse1]['tilt']['x']
                                sse1_norm = abs(sse1_tx)/sse1_tx
                                sse2_tx   = data['corrections'][sse2]['tilt']['x']
                                sse2_norm = abs(sse2_tx)/sse2_tx
                                if sse1_norm != sse2_norm:
                                    self.log.debug(f'Re-correcting {sse2} tilt.')
                                    data['corrections'][sse2]['tilt']['x'] *= -1.

        self.log.notice(f'Found corrections {data["corrections"]}\n')
        return data


    def alpha_on_beta_correction(self, df: pd.DataFrame, bin: str, wdir: Path, qlayer: str, rlayer: str,
                                 case: Case, extras: pd.DataFrame ) -> List[Dict]:
        """
        """
        # Report data
        stats = self.make_mode_stats(df, wdir).reset_index()
        clms = ['measure', 'layer', 'sse', bin]
        stats = stats[(stats['layer'] == rlayer)][clms]
        extras = extras[(extras['layer'] == rlayer)]
        for layer in sorted(df.layer.unique()):
            ofile = 'geometric_distributions_layer{}'.format(layer)
            TButil.plot_geometric_distributions(self.log, df[df['layer'] == layer], Path(wdir).joinpath(ofile))

        data = {}
        preref = {'angles_layer': 0, 'angles_side': 0, 'angles_floor': 0,
                  'points_side': 0, 'points_side': 0, 'points_layer': 0}
        for sse in [x for x in stats.sse.unique() if x.startswith(qlayer)]:
            # pre-calculations here
            ddf = stats[(stats['sse'] == sse)]
            ddx = extras[(extras['sse'] == sse)]
            #preref['angles_layer'] = -ddx['angles_layer'].values[0]
            #preref['angles_side'] = -ddx['angles_side'].values[0]
            #preref['angles_floor'] = -ddx['angles_floor'].values[0]

            d_tilt, d_trans = {}, {}
            for key, val in self.correctives.items():
                self.log.debug(key, val)
                if key == 'rotate':
                    for v in val:
                        if v == 'x':
                            preref['angles_layer'] = -ddx['angles_layer'].values[0]
                            d_tilt[v] = ddf[ddf['measure'] == 'angles_layer'][bin].values[0] + preref['angles_layer']
                        elif v == 'y':
                            preref['angles_floor'] = -ddx['angles_floor'].values[0]
                            d_tilt[v] = ddf[ddf['measure'] == 'angles_floor'][bin].values[0] + preref['angles_floor']
                        elif v == 'z':
                            preref['angles_side'] = -ddx['angles_side'].values[0]
                            d_tilt[v] = ddf[ddf['measure'] == 'angles_side'][bin].values[0] + preref['angles_side']
                        else:
                            self.log.info(f'Unknown rotation for {v}, not applying anything.')
                            #raise NodeOptionsError(f'Unknown rotation for {v}')
                if key == 'translate':
                    for v in val:
                        if v == 'x':
                            preref['points_side'] = -ddx['points_side'].values[0]
                            d_trans[v] = ddf[ddf['measure'] == 'points_side'][bin].values[0] + preref['points_side']
                        elif v == 'y':
                            preref['points_floor'] = -ddx['points_floor'].values[0]
                            d_trans[v] = ddf[ddf['measure'] == 'points_floor'][bin].values[0] + preref['points_floor']
                        elif v == 'z':
                            self.log.debug(case['configuration.defaults.distance.ab'])
                            self.log.debug(ddf[ddf['measure'] == 'points_layer'][bin].values[0])
                            pc = case['configuration.defaults.distance.ab'] - ddf[ddf['measure'] == 'points_layer'][bin].values[0]
                            if ascii_uppercase.index(qlayer) < ascii_uppercase.index(rlayer):
                                pc = pc * -1
                            d_trans[v] = pc
                        else:
                            self.log.info(f'Unknown translation for {v}, not applying anything.')
                            #raise NodeOptionsError(f'Unknown translation for {v}')

            data.setdefault(sse, {}).setdefault('tilt', d_tilt)
            data.setdefault(sse, {}).setdefault('coordinates', d_trans)

            #data.setdefault(sse, {}).setdefault('tilt', {'x': ddf[ddf['measure'] == 'angles_layer'][bin].values[0] + preref['angles_layer'],
            #                                             'z': ddf[ddf['measure'] == 'angles_side'][bin].values[0] + preref['angles_side']})
            #data.setdefault(sse, {}).setdefault('tilt', {'x': ddf[ddf['measure'] == 'angles_layer'][bin].values[0] + preref['angles_layer'],})
            #data.setdefault(sse, {}).setdefault('tilt', {'z': ddf[ddf['measure'] == 'angles_side'][bin].values[0] + preref['angles_side'],})

            #pc = ddf[ddf['measure'] == 'points_layer'][bin].values[0] - case['configuration.defaults.distance.ab']
            #if ascii_uppercase.index(qlayer) < ascii_uppercase.index(rlayer):
            #    pc = pc * -1
            #preref['points_floor'] = -ddx['points_floor'].values[0]

            #data.setdefault(sse, {}).setdefault('coordinates', {'y': ddf[ddf['measure'] == 'points_floor'][bin].values[0] + preref['points_floor'],
            #                                                    'z': pc})
            #data.setdefault(sse, {}).setdefault('coordinates', {'y': ddf[ddf['measure'] == 'points_floor'][bin].values[0] + preref['points_floor'],})
            #data.setdefault(sse, {}).setdefault('coordinates', {'z': pc,})
        return data, preref


    def beta_on_beta_correction(self, df: pd.DataFrame, bin: str, wdir: Path, qlayer: str, rlayer: str,
                                 case: Case, extras: pd.DataFrame, rules: list ) -> List[Dict]:
        """
        """
        # Report data
        stats = self.make_mode_stats(df, wdir).reset_index()
        clms = ['measure', 'layer', 'sse', bin]
        stats = stats[(stats['layer'] == qlayer)][clms]
        extras = extras[(extras['layer'] == qlayer)]
        for layer in sorted(df.layer.unique()):
            ofile = 'geometric_distributions_layer{}'.format(layer)
            TButil.plot_geometric_distributions(self.log, df[df['layer'] == layer], Path(wdir).joinpath(ofile))

        data = {}
        preref = {'angles_layer': 0, 'angles_side': 0, 'angles_floor': 0,
                  'points_side': 0, 'points_side': 0, 'points_layer': 0}

        angle_side      = stats[(stats['measure'] == 'angles_side')][bin]
        angle_side_sign = (angle_side.abs() / angle_side).values
        angle_side_mean = angle_side.abs().mean() * (np.abs(angle_side_sign.sum()) / angle_side_sign.sum())
        for sse in [x for x in stats.sse.unique() if x.startswith(qlayer)]:
            for rule in rules:
                if rule[0] == sse:
                    if rule[-1] == True:
                        flip = -1.
                    else:
                        flip = 1.

            ddf = stats[(stats['sse'] == sse)]
            ddx = extras[(extras['sse'] == sse)]
            #preref['angles_layer'] = -ddx['angles_layer'].values[0]
            #preref['angles_side'] = -ddx['angles_side'].values[0]

            d_tilt, d_trans = {}, {}
            for key, val in self.correctives.items():
                if key == 'rotate':
                    for v in val:
                        if v == 'x':
                            preref['angles_layer'] = -ddx['angles_layer'].values[0]
                            d_tilt[v] = ddf[ddf['measure'] == 'angles_layer'][bin].values[0] + preref['angles_layer']
                        elif v == 'y':
                            preref['angles_floor'] = -ddx['angles_floor'].values[0]
                            d_tilt[v] = ddf[ddf['measure'] == 'angles_floor'][bin].values[0] + preref['angles_floor']
                        elif v == 'z':
                            #preref['angles_side'] = -ddx['angles_side'].values[0]
                            d_tilt[v] = angle_side_mean * flip
                        else:
                            self.log.info(f'Unknown rotation for {v}, not applying anything.')
                            #raise NodeOptionsError(f'Unknown rotation for {v}')
                if key == 'translate':
                    for v in val:
                        if v == 'x':
                            preref['points_side'] = -ddx['points_side'].values[0]
                            d_trans[v] = ddf[ddf['measure'] == 'points_side'][bin].values[0] + preref['points_side']
                        elif v == 'y':
                            preref['points_floor'] = -ddx['points_floor'].values[0]
                            d_trans[v] = ddf[ddf['measure'] == 'points_floor'][bin].values[0] + preref['points_floor']
                        elif v == 'z':
                            #pc = case['configuration.defaults.distance.bb_stack'] + ddf[ddf['measure'] == 'points_layer'][bin].values[0]
                            pc = ddf[ddf['measure'] == 'points_layer'][bin].values[0]
                            if ascii_uppercase.index(qlayer) < ascii_uppercase.index(rlayer):
                                pc = pc * -1.
                            d_trans[v] = pc
                        else:
                            self.log.info(f'Unknown translation for {v}, not applying anything.')
                            raise NodeOptionsError(f'Unknown translation for {v}')

            data.setdefault(sse, {}).setdefault('tilt', d_tilt)
            data.setdefault(sse, {}).setdefault('coordinates', d_trans)

            #data.setdefault(sse, {}).setdefault('tilt', {'x': ddf[ddf['measure'] == 'angles_layer'][bin].values[0] + preref['angles_layer'],
            #                                             'z': angle_side_mean * flip})
            #data.setdefault(sse, {}).setdefault('tilt', {'x': ddf[ddf['measure'] == 'angles_layer'][bin].values[0] + preref['angles_layer'],})
            #data.setdefault(sse, {}).setdefault('tilt', {'z': angle_side_mean * flip,})

            #pc = ddf[ddf['measure'] == 'points_layer'][bin].values[0] #- case['configuration.defaults.distance.bb_stack']
            #if ascii_uppercase.index(qlayer) < ascii_uppercase.index(rlayer):
            #   pc = pc * -1
            #preref['points_floor'] = -ddx['points_floor'].values[0]

            #data.setdefault(sse, {}).setdefault('coordinates', {'y': ddf[ddf['measure'] == 'points_floor'][bin].values[0] + preref['points_floor'],
            #                                                    'z': pc })
            #data.setdefault(sse, {}).setdefault('coordinates', {'y': ddf[ddf['measure'] == 'points_floor'][bin].values[0] + preref['points_floor'],})
            #data.setdefault(sse, {}).setdefault('coordinates', {'z': pc,})
        return data, preref


    def nth_layer_correction( self, df: pd.DataFrame, case: Case, bin: str, wdir: Path, qlayer: str, rlayer: str, rules: list) -> Dict:
        """
        """
        def create_network(df, bin, wdir, qlayer, rlayer, corr_type):
            def reshape(df):
                def reshape(g):
                    g = g.drop(columns=['pdb', 'chain']).T
                    g = g.rename(columns=g.loc['sse'])
                    g = g.reindex(g.index.drop('sse'))
                    g = g.assign(bin=g.loc['bin'].values[0])
                    g = g.assign(N0X=0)
                    g = g.drop(index='bin')
                    return g
                df = df.copy()[['pdb', 'chain', corr_type, 'sse', 'bin']]
                bins = list(range(-100, 105, 5))
                labels = list(np.arange(-97.5, 100, 5))
                df = df.assign(anglebin=pd.cut(df[corr_type], bins=bins, labels=labels))
                df = df.drop_duplicates(["pdb", "chain", "sse", "bin"])
                df = df.drop(columns=[corr_type]).groupby(['pdb', 'chain']).apply(reshape)
                df.index = list(range(df.shape[0]))
                return df.drop(columns=['bin'])
            ddf = reshape(df[(df['bin'] == bin) & (df['layer'] == qlayer)])
            sses = sorted(list(ddf.columns))
            netwk = []
            for i in range(0, len(sses) - 1):
                idx = [sses[i], sses[i + 1]]
                tmp = pd.DataFrame(ddf.groupby(idx).size()).rename(columns={0: 'count'}).reset_index()
                netwk.append(tmp.reset_index())
                netwk[-1][idx[0]] = netwk[-1][idx[0]].apply(lambda v: "_".join([idx[0], str(v)]))
                netwk[-1][idx[1]] = netwk[-1][idx[1]].apply(lambda v: "_".join([idx[1], str(v)]))
                netwk[-1] = nx.from_pandas_edgelist(netwk[-1], idx[0], idx[1], ['count'], create_using=nx.DiGraph)

            posk = {}
            try:
                network = nx.compose_all(netwk)
                for n in network.nodes:
                    dd = n.split('_')
                    posk.setdefault(n, (int(dd[0][1]), float(dd[1])))
            except Exception:
                network = nx.Graph()

            # Image summary
            TButil.plot_angle_network(self.log, network, posk, sses, Path(wdir).joinpath(f'network_{corr_type}_{bin}'))
            return network

        # Report data
        self.make_mode_stats(df, wdir)
        for layer in sorted(df.layer.unique()):
            ofile = 'geometric_distributions_layer{}'.format(layer)
            TButil.plot_geometric_distributions(self.log, df[df['layer'] == layer], Path(wdir).joinpath(ofile))

        # Rotation networks
        keys = df['sse'].drop_duplicates().tolist()
        d_tilts = {k: {} for k in keys}
        d_trans = {k: {} for k in keys}
        for key, val in self.correctives.items():
            if key == 'rotate':
                for v in val:
                    if v == 'x':
                        networkrx = create_network(df, bin, wdir, qlayer, rlayer, 'angles_layer')
                        for e in nx.dag_longest_path(networkrx, 'count', default_weight=0):
                            e = e.split('_')
                            if e[0] == 'N0X':
                                continue
                            elif e[0].startswith(qlayer):
                                #d_tilt[v] = float(e[1])
                                d_tilts[e[0]].update({v: float(e[1])})
                            else:
                                continue
                    elif v == 'y':
                        networkry = create_network(df, bin, wdir, qlayer, rlayer, 'angles_floor')
                        for e in nx.dag_longest_path(networkry, 'count', default_weight=0):
                            #d_tilt = {}
                            e = e.split('_')
                            if e[0] == 'N0X':
                                continue
                            elif e[0].startswith(qlayer):
                                #d_tilt[v] = float(e[1])
                                d_tilts[e[0]].update({v: float(e[1])})
                            else:
                                continue
                    elif v == 'z':
                        networkrz = create_network(df, bin, wdir, qlayer, rlayer, 'angles_side')
                        tiltsz, tiltz_sses = [], []
                        for e in nx.dag_longest_path(networkrz, 'count', default_weight=0):
                            #d_tilt = {}
                            e = e.split('_')
                            if e[0] == 'N0X':
                                continue
                            elif e[0].startswith(qlayer):
                                #d_tilt[v] = float(e[1])
                                #d_tilts[e[0]].update({v: float(e[1])})
                                tiltsz.append(float(e[1]))
                                tiltz_sses.append(e[0])
                            else:
                                continue
                        tiltsz_mean = abs(np.mean(tiltsz))
                        tiltsz_sum  = np.sum(tiltsz)
                        orienter = tiltsz_sum / abs(tiltsz_sum)
                        for ii in range(len(tiltz_sses)):
                            for rule in rules:
                                if rule[0] == tiltz_sses[ii]:
                                    if rule[-1] == True:
                                        if orienter == -1.:
                                            flip = 1.
                                        else:
                                            flip = -1.
                                    else:
                                        if orienter == -1.:
                                            flip = -1.
                                        else:
                                            flip = 1.
                            d_tilts[tiltz_sses[ii]].update({v: float(tiltsz_mean * flip)})
                    else:
                        self.log.info(f'Unknown rotation for {v}, not applying anything.')
                        #raise NodeOptionsError(f'Unknown rotation for {v}')
            if key == 'translate':
                for v in val:
                    if v == 'x':
                        networktx = create_network(df, bin, wdir, qlayer, rlayer, 'points_side')
                        for e in nx.dag_longest_path(networktx, 'count', default_weight=0):
                            #d_tran = {}
                            e = e.split('_')
                            if e[0] == 'N0X':
                                continue
                            elif e[0].startswith(qlayer):
                                #d_tran[v] = float(e[1])
                                d_trans[e[0]].update({v: float(e[1])})
                            else:
                                continue
                    elif v == 'y':
                        networkty = create_network(df, bin, wdir, qlayer, rlayer, 'points_floor')
                        for e in nx.dag_longest_path(networkty, 'count', default_weight=0):
                            #d_tran = {}
                            e = e.split('_')
                            if e[0] == 'N0X':
                                continue
                            elif e[0].startswith(qlayer):
                                #d_tran[v] = float(e[1])
                                d_trans[e[0]].update({v: float(e[1])})
                            else:
                                continue
                    elif v == 'z':
                        networktz = create_network(df, bin, wdir, qlayer, rlayer, 'points_layer')
                        for e in nx.dag_longest_path(networktz, 'count', default_weight=0):
                            #d_tran = {}
                            e = e.split('_')
                            if e[0] == 'N0X':
                                continue
                            elif e[0].startswith(qlayer):
                                pc = float(e[1]) #- case['configuration.defaults.distance.ab']
                                if ascii_uppercase.index(qlayer) < ascii_uppercase.index(rlayer):
                                    pc = pc * -1.
                                #d_tran[v] = pc
                                d_trans[e[0]].update({v: float(e[1])})
                            else:
                                continue
                    else:
                        self.log.info(f'Unknown translation for {v}, not applying anything.')
                        #raise NodeOptionsError(f'Unknown translation for {v}')
        d_trans = {k: v for k,v in d_trans.items() if v != {}}
        d_tilts = {k: v for k,v in d_tilts.items() if v != {}}
        data = {}
        for k,v in d_trans.items():
            data.setdefault(k, {}).setdefault('coordinates', v)
        for k,v in d_tilts.items():
            data.setdefault(k, {}).setdefault('tilt', v)
        return data


    def first_layer_correction( self, df: pd.DataFrame, bin: str, wdir: Path ) -> Dict:
        """
        """
        # Report data
        self.make_mode_stats(df, wdir)
        for layer in sorted(df.layer.unique()):
            ofile = 'geometric_distributions_layer{}'.format(layer)
            TButil.plot_geometric_distributions(self.log, df[df['layer'] == layer], Path(wdir).joinpath(ofile))

        # Make network
        def reshape(df):
            def reshape(g):
                g = g.drop(columns=['pdb', 'chain']).T
                g = g.rename(columns=g.loc['sse'])
                g = g.reindex(g.index.drop('sse'))
                g = g.assign(bin=g.loc['bin'].values[0])
                g = g.assign(A0E=0)
                g = g.drop(index='bin')
                return g

            df = df.copy()[['pdb', 'chain', 'angles_layer', 'sse', 'bin']]
            bins = list(range(-100, 105, 5))
            labels = list(np.arange(-97.5, 100, 5))
            df = df.assign(anglebin=pd.cut(df['angles_layer'], bins=bins, labels=labels))
            df = df.drop_duplicates(["pdb", "chain", "sse", "bin"])
            df = df.drop(columns=['angles_layer']).groupby(['pdb', 'chain']).apply(reshape)
            df.index = list(range(df.shape[0]))
            return df.drop(columns=['bin'])
        ddf = reshape(df[df['bin'] == bin])
        sses = sorted(list(ddf.columns))
        netwk = []
        for i in range(0, len(sses) - 1):
            idx = [sses[i], sses[i + 1]]
            tmp = pd.DataFrame(ddf.groupby(idx).size()).rename(columns={0: 'count'}).reset_index()
            netwk.append(tmp.reset_index())
            netwk[-1][idx[0]] = netwk[-1][idx[0]].apply(lambda v: "_".join([idx[0], str(v)]))
            netwk[-1][idx[1]] = netwk[-1][idx[1]].apply(lambda v: "_".join([idx[1], str(v)]))
            netwk[-1] = nx.from_pandas_edgelist(netwk[-1], idx[0], idx[1], ['count'], create_using=nx.DiGraph)

        posk = {}
        try:
            network = nx.compose_all(netwk)
            for n in network.nodes:
                dd = n.split('_')
                posk.setdefault(n, (int(dd[0][1]), float(dd[1])))
        except Exception:
            network = nx.Graph()

        # Image summary
        TButil.plot_angle_network(self.log, network, posk, sses, Path(wdir).joinpath('network_{}'.format(bin)))

        # Return corrections
        data = {}
        for x in nx.dag_longest_path(network, 'count', default_weight=0):
            x = x.split('_')
            if x[0] == 'A0E':
                continue
            else:
                data.setdefault(x[0], {}).setdefault('tilt', {}).setdefault('x', float(x[1]))
        return data


    def first_layer_alpha_correction( self, df: pd.DataFrame, bin: str, wdir: Path ) -> Dict:
        """
        """
        # Report data
        self.make_mode_stats(df, wdir)
        for layer in sorted(df.layer.unique()):
            ofile = 'geometric_distributions_layer{}'.format(layer)
            TButil.plot_geometric_distributions(self.log, df[df['layer'] == layer], Path(wdir).joinpath(ofile))

        # Make network
        def reshape(df):
            def reshape(g):
                g = g.drop(columns=['pdb', 'chain']).T
                g = g.rename(columns=g.loc['sse'])
                g = g.reindex(g.index.drop('sse'))
                g = g.assign(bin=g.loc['bin'].values[0])
                g = g.assign(A0E=0)
                g = g.drop(index='bin')
                return g

            df = df.copy()[['pdb', 'chain', 'angles_layer', 'sse', 'bin']]
            bins = list(range(-100, 105, 5))
            labels = list(np.arange(-97.5, 100, 5))
            df = df.assign(anglebin=pd.cut(df['angles_layer'], bins=bins, labels=labels))
            df = df.drop(columns=['angles_layer']).groupby(['pdb', 'chain']).apply(reshape)
            df.index = list(range(df.shape[0]))
            return df.drop(columns=['bin'])
        ddf = reshape(df[df['bin'] == bin])
        sses = sorted(list(ddf.columns))
        netwk = []
        for i in range(0, len(sses) - 1):
            idx = [sses[i], sses[i + 1]]
            tmp = pd.DataFrame(ddf.groupby(idx).size()).rename(columns={0: 'count'}).reset_index()
            netwk.append(tmp.reset_index())
            netwk[-1][idx[0]] = netwk[-1][idx[0]].apply(lambda v: "_".join([idx[0], str(v)]))
            netwk[-1][idx[1]] = netwk[-1][idx[1]].apply(lambda v: "_".join([idx[1], str(v)]))
            netwk[-1] = nx.from_pandas_edgelist(netwk[-1], idx[0], idx[1], ['count'], create_using=nx.DiGraph)

        posk = {}
        try:
            network = nx.compose_all(netwk)
            for n in network.nodes:
                dd = n.split('_')
                posk.setdefault(n, (int(dd[0][1]), float(dd[1])))
        except Exception:
            network = nx.Graph()

        # Image summary
        TButil.plot_angle_network(self.log, network, posk, sses, Path(wdir).joinpath('network_{}'.format(bin)))

        # Return corrections
        data = {}
        for x in nx.dag_longest_path(network, 'count', default_weight=0):
            x = x.split('_')
            if x[0] == 'A0H':
                continue
            else:
                data.setdefault(x[0], {}).setdefault('tilt', {}).setdefault('x', float(x[1]))
        return data


    def make_mode_stats( self, df: pd.DataFrame, wdir: Path ) -> pd.DataFrame:
        """
        """
        data = {'close': [], 'mid': [], 'far': [], 'extreme': []}
        topi, midi, boti = [], [], []

        for angle in [x for x in df.columns if x.startswith('angles_') or x.startswith('points_')]:
            for part in ['close', 'mid', 'far', 'extreme']:
                dfp = df[df['bin'] == part]
                for lay in df.layer.unique():
                    for sse in df.sse.unique():
                        dfs = dfp[(dfp['sse'] == sse) & (dfp['layer'] == lay)]
                        if not dfs[angle].empty:
                            try:
                                kde = sc.stats.gaussian_kde(dfs[angle])
                                x = np.linspace(dfs[angle].min(), dfs[angle].max(), 200)
                                kde = kde(x)
                                mode = x[np.argsort(kde)[-1]]
                            except ValueError:
                                mode = dfs[angle].values[-1]
                        else:
                            mode = 0.
                        data[part].append(mode)
                        if part == 'close':
                            topi.append(angle)
                            midi.append(lay)
                            boti.append(sse)
        stats = pd.DataFrame(data, index=pd.MultiIndex.from_tuples(list(zip(*[topi, midi, boti])),
                                                                   names=['measure', 'layer', 'sse']))
        stats.reset_index().to_csv(wdir.joinpath('mode_stats.csv'), index=False)
        self.log.notice('Mode stats stored at {}'.format(wdir.joinpath('mode_stats.csv')))
        return stats


    def with_slurm( self, cmd: List[str],
                    current_case_file: Path,
                    current_sse: str,
                    unimaster: Path,
                    imaster: Path,
                    unidata: Path ):
        """
        """
        # Make bashfile
        bashcont = []
        createbash = 'python {0} -case {1} -master {2} -present {3} -out {4}'
        parts = math.ceil(len(cmd) / TBcore.get_option('slurm', 'array'))

        wwd = unimaster.parent.parent
        cwd = Path().cwd()
        os.chdir(str(wwd))

        for i, com in enumerate(cmd):
            cmd[i][2] = str(Path(com[2]).relative_to(wwd))
            cmd[i][-1] = str(Path(com[-1]).relative_to(wwd))

        for ii, cp in enumerate(cmd):
            cmd[ii][-1] = cp[-1] + '_${SLURM_ARRAY_TASK_ID}'
        for j, i in enumerate(range(0, len(cmd), parts)):
            sumfile = unimaster.parent.joinpath('_${SLURM_ARRAY_TASK_ID}.master').relative_to(wwd)
            datfile = unimaster.parent.joinpath('_${SLURM_ARRAY_TASK_ID}.geo').relative_to(wwd)
            bashcont.append('if (( ${{SLURM_ARRAY_TASK_ID}} == {} )); then'.format(j + 1))
            bashcont.extend([' '.join(x) for x in cmd[i:i + parts]])
            bashcont.append('cat {0} > {1}'.format(Path(cmd[-1][-1]).parent.joinpath('*_${SLURM_ARRAY_TASK_ID}'), sumfile))
            bashcont.append(createbash.format(imaster, current_case_file.relative_to(wwd),
                                              sumfile, current_sse, datfile))
            bashcont.append('fi')
        with unimaster.parent.joinpath('submit.sh').relative_to(wwd).open('w') as fd:
            fd.write(TButil.slurm_header())
            fd.write(TButil.slurm_pyenv())
            # due to slurms goddam size limits, we need to keep the submitting
            # file short and put all the lines in an other file.
            #fd.write(TButil.slurm_exec())
            fd.write('\n' + 'bash {} ${{SLURM_ARRAY_TASK_ID}}'.format(
                unimaster.parent.joinpath('exec.sh').relative_to(wwd)) + '\n')
        with unimaster.parent.joinpath('exec.sh').relative_to(wwd).open('w') as fd:
            fd.write(TButil.bash_exec_header())
            fd.write('\n'.join(bashcont))

        TButil.submit_slurm(self.log, unimaster.parent.joinpath('submit.sh').relative_to(wwd))
        self.log.notice(f'Creating geometric coordinate file {unidata}')
        allCSV = [str(x) for x in unimaster.parent.relative_to(wwd).glob('_*.geo.csv')]
        pd.concat([pd.read_csv(x) for x in allCSV]).to_csv(unidata.relative_to(wwd), index=False)
        self.log.notice(f'Creating MASTER search file {unimaster}')
        with unimaster.relative_to(wwd).open('w') as fd:
            for x in unimaster.parent.glob('_*.master'):
                with x.relative_to(wwd).open() as fi:
                    fd.write(''.join(fi.readlines()))
        os.chdir(str(cwd))


    def no_slurm( self, cmd: List[str],
                  current_case_file: Path,
                  current_sse: str,
                  unimaster: Path,
                  imaster: Path,
                  unidata: Path):
        """
        """
        # Search on MASTER
        result = []
        for com in cmd:
            self.log.notice(f'EXECUTE: {" ".join([str(x) for x in com])}')
            run(com, stdout=DEVNULL)
            outf = Path(com[-1])
            if outf.is_file():
                result.append(str(outf))

        try:
            result.insert(0, 'cat')
            self.log.notice(f'Unify matches at {unimaster}')
            with unimaster.open('w') as fd:
                run(result, stdout=fd)
        except:
            with unimaster.open('w') as fd:
                for ms in glob.iglob(os.path.dirname(result[-1]) + '/*.master'):
                    with open(ms, 'r') as f:
                        lines = f.readlines()
                    for line in lines:
                        fd.write(line)

        if self.subsampling: # random subsampling due to local (no_slurm)
            self.log.debug(f'Sub-sampling {self.subsampling} from full master matches.')
            with unimaster.open('r') as fd:
                lines = fd.readlines()
            lines = np.random.choice(lines, self.subsampling, replace=False)
            unimastersub = str(unimaster).replace('match.master', 'matchsub.master')
            with open(unimastersub, 'w') as f:
                for line in lines:
                    f.write(line)
            # Analyze
            createbash = 'python {0} -case {1} -master {2} -present {3} -out {4}'
            cmd = shlex.split(createbash.format(imaster, current_case_file, unimastersub, current_sse, unidata))
            self.log.notice(f'EXECUTE: {" ".join([str(x) for x in cmd])}')
            run(cmd, stdout=DEVNULL)

        else:
            createbash = 'python {0} -case {1} -master {2} -present {3} -out {4}'
            cmd = shlex.split(createbash.format(imaster, current_case_file, unimaster, current_sse, unidata))
            self.log.notice(f'EXECUTE: {" ".join([str(x) for x in cmd])}')
            run(cmd, stdout=DEVNULL)
