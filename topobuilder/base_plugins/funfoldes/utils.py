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
import math
import gzip
from pathlib import Path
from logbook import Logger
from typing import Union, Dict, Tuple, List, Optional
import sys
import textwrap
from subprocess import run, DEVNULL

# External Libraries
import pandas as pd
from rstoolbox.io import open_rosetta_file
try:
    from SBI.structure import PDB, Frame3D
    import SBI.core as SBIcr
except ImportError:
    class Frame3D():
        pass

# This Library
from topobuilder.case import Case
import topobuilder.core as TBcore
import topobuilder.utils as TButil
from topobuilder.utils import build_pdb_object
from topobuilder.workflow import TBplugins


def folder_structure( case: Case ) -> Dict:
    """Create the folder structure of the plugin.
    """
    # Generate the folder tree for a single connectivity.
    folders = case.connectivities_paths[0].joinpath('funfoldes')
    folders.mkdir(parents=True, exist_ok=True)
    outdir = folders.joinpath('outputs')
    outdir.mkdir(parents=True, exist_ok=True)
    pdb_file = folders.joinpath('template_sketch.pdb')
    ffd_fold_file = folders.joinpath('funfoldes_fold.xml')
    ffd_design_file = folders.joinpath('funfoldes_design.xml')
    wts0_file = folders.joinpath('score0.wts')
    wts0_patch_file = folders.joinpath('score0.wts_patch')
    wts1_file = folders.joinpath('score1.wts')
    wts1_patch_file = folders.joinpath('score1.wts_patch')
    wts2_file = folders.joinpath('score2.wts')
    wts2_patch_file = folders.joinpath('score2.wts_patch')
    wts3_file = folders.joinpath('score3.wts')
    wts3_patch_file = folders.joinpath('score3.wts_patch')
    wts5_file = folders.joinpath('score5.wts')
    wts5_patch_file = folders.joinpath('score5.wts_patch')
    checkpoint = folders.joinpath('checkpoint.json')

    return {'main': folders,                # Main plugin folder
            'outdir': outdir,               # Folder to produce the Rosetta outputs
            'pdb': pdb_file,                # Path to the template PDB file
            'foldRS': ffd_fold_file,        # Path to the Fold script
            'designRS': ffd_design_file,    # Path to the Design script
            'wts0F': wts0_file,
            'wts0_patchF': wts0_patch_file,
            'wts1F': wts1_file,
            'wts1_patchF': wts1_patch_file,
            'wts2F': wts2_file,
            'wts2_patchF': wts2_patch_file,
            'wts3F': wts3_file,
            'wts3_patchF': wts3_patch_file,
            'wts5F': wts5_file,
            'wts5_patchF': wts5_patch_file,
            'checkpoint': checkpoint        # Path to the Checkpoint file (.json)
            }


def build_template_sketch( log: Logger, case: Case, full_file: Union[Path, str] ):
    """Generate the PDB file used as template to fold.
    """
    bindersfile = None
    pdb_file = os.path.dirname( str(full_file) ) + f'/motif.pdb'
    columns = ['auth_comp_id', 'auth_atom_id', 'auth_asym_id', 'auth_seq_id', 'Cartn_x', 'Cartn_y', 'Cartn_z']
    # List setup for saving relevant info if necessary
    #Spots = namedtuple('Spots', ['motifs', 'binders', 'hotpsots', 'coldspots', 'identifier'])
    res_attach, res_hotspots, res_coldspots, binder_chains, m_identifiers = [], [], [], [], []

    # 1.1 Make sure the loop_length are specified.
    loop_lengths = case['metadata.loop_lengths']
    looplist = ', '.join([str(x) for x in loop_lengths])
    log.info(f'Gathering a total of {len(loop_lengths)} loops of length: {looplist}\n')
    log.info(f'To apply over {len(case)} secondary structures: {case.connectivities_str[0]}\n')

    # Get the structure.
    #if not case.data['metadata']['motif_picker']: # no motif
    if not 'motif_picker' in case.data['metadata']:
        node = getattr(TBplugins.source.load_plugin('builder'), 'builder', None)(connectivity=True, pick_aa='V', tag=0)
    else: # motif
        node = getattr(TBplugins.source.load_plugin('builder'), 'builder', None)(connectivity=True, motif=True, pick_aa='V', tag=0)
    case = Case(node.single_execute(case.data))

    # Get the binder.
    #if not case.data['metadata']['motif_picker']: # no motif
    if not 'motif_picker' in case.data['metadata']:
        node = getattr(TBplugins.source.load_plugin('builder'), 'builder', None)(connectivity=True, pick_aa='V', tag=0)
    else: # motif and hotspots
        node = getattr(TBplugins.source.load_plugin('builder'), 'builder', None)(connectivity=True, motif=True, pick_aa='V', tag=0)
    case = Case(node.single_execute(case.data))

    sse_list = case.ordered_structures
    pdb, _ = build_pdb_object(log, sse_list, loop_lengths)
    log.notice(f'Writing structure {pdb_file}')
    pdb.write(str(pdb_file), format='pdb', clean=True, force=TBcore.get_option('system', 'overwrite'))

    #if not case.data['metadata']['binder']:
    if not 'binder' in case.data['metadata']:
        log.notice(f'Writing structure (no binder) {full_file}')
        pdb.write(str(full_file), format='pdb', clean=True, force=TBcore.get_option('system', 'overwrite'))

    #if case.data['metadata']['motif_picker']:
    if 'motif_picker' in case.data['metadata']:
        for mdata in case.data['metadata']['motif_picker']:
            motif, binder, hotspots, attach, selection, identifier = mdata['motifs']

            attach = pdb[pdb.sse_id.isin(attach)][['auth_seq_id', 'auth_asym_id']].drop_duplicates().values
            attach = [f'{s[0]}{s[1]}' for s in attach]
            hotspots  = pdb[pdb.pdb_num.isin(hotspots['auth_seq_id'].drop_duplicates().values
                            )][['auth_seq_id', 'auth_asym_id']].drop_duplicates().values
            hotspots  = [f'{s[0]}{s[1]}' for s in hotspots]
            coldspots = [s for s in attach if s not in hotspots]
            log.debug(f'Adding motif of id {identifier}: {attach}')
            log.debug(f'Current hotspots  : {hotspots}')
            log.debug(f'Current coldspots : {coldspots}')

            res_attach.extend(attach)
            res_hotspots.extend(hotspots)
            res_coldspots.extend(coldspots)
            m_identifiers.extend(identifier)

    #if case.data['metadata']['binder']: # Binder
    if 'binder' in case.data['metadata']:
        #full_structure = [pdb,]
        binders = []
        for key in case.data['metadata']['binder']:
            binder = case.data['metadata']['binder'][key]
            binders.append(binder)
            binderfile = os.path.dirname( str(pdb_file) ) + f'/binder_{key}.pdb'
            binder_chains.extend( binder['auth_asym_id'].drop_duplicates().tolist() )
            #full_structure.append(binder)
            log.debug(f'Adding binder chains: {binder_chains}')

        bindersfile = os.path.dirname( str(pdb_file) ) + f'/binders.pdb'
        log.notice(f'Writing structure {bindersfile}')
        binders = PDB(pd.concat(binders, sort=False))
        binders.write(bindersfile, format='pdb', clean=True, force=TBcore.get_option('system', 'overwrite'))

        full_structure = PDB( pd.concat([pdb[columns], binders[columns]], sort=False) )
        log.notice(f'Writing structure {full_file}')
        full_structure.write(str(full_file), format='pdb', clean=True, force=TBcore.get_option('system', 'overwrite'))
    #else:
    #    pdb.write(str(pdb_file), format='pdb', clean=True, force=TBcore.get_option('system', 'overwrite'))

    # Push back hotspots, motif seq num and binder chain
    if binder_chains == []: binder_chains = None
    if res_attach == []:    res_attach = None
    if res_hotspots == []:  res_hotspots = None

    return binder_chains, res_attach, res_hotspots, res_coldspots, m_identifiers, bindersfile


def make_scripts( log: Logger,
                  case: Case,
                  wpaths: Dict,
                  data: Dict,
                  natbias: float = 2.5,
                  layer_design: bool = True,
                  binder: Optional[List] = None,
                  motif: Optional[List] = None,
                  hotspots: Optional[List] = None,
                  identifiers: Optional[List] = None,
                  binderfile: Optional[str] = None,
                  profile: Optional[bool] = False,
                  ) -> Tuple[str, str]:
    """Create the folding and design scripts.
    """
    fld = TButil.rosettascript(TButil.funfoldes(case, motif, binder, hotspots))
    dsg = TButil.rosettascript(TButil.constraint_design(case, natbias, layer_design,
                                                        motif, binder, hotspots, profile))
    wts = TButil.get_weight_patches()

    if TBcore.get_option('system', 'jupyter'):
        ifold = os.getenv('TB_FUNFOLDES_FOLD_FILE', None)
        idsgn = os.getenv('TB_FUNFOLDES_DSGN_FILE', None)

        if ifold is None:
            print('-' * 80)
            print('\n\nExpected FOLDING script will be:\n')
            print(fld)
            print(textwrap.dedent("""\n\n If a different script wants to be provided, do so by
    assigning a RosettaScript to os.environ['TB_FUNFOLDES_FOLD_FILE'] or assign '' to
    it if you are ok with the default script."""))
        elif ifold:
            ifold = Path(ifold)
            if not ifold.is_file():
                raise IOError(f'Unknown file {ifold}')
            fld = ''.join(list(ifold.open().readlines()))

        if idsgn is None:
            print('-' * 80)
            print('\n\nExpected DESIGN script will be:\n')
            print(dsg)
            print(textwrap.dedent("""\n\n If a different script wants to be provided, do so by
    assigning a RosettaScript to os.environ['TB_FUNFOLDES_DSGN_FILE'] or assign '' to
    it if you are ok with the default script."""))
        elif idsgn:
            idsgn = Path(idsgn)
            if not idsgn.is_file():
                raise IOError(f'Unknown file {idsgn}')
            dsg = ''.join(list(idsgn.open().readlines()))

        if ifold is None or idsgn is None:
            TButil.exit()

    log.info(f'Writing the folding RosettaScript file: {wpaths["foldRS"]}\n')
    with wpaths['foldRS'].open('w') as fd:
        fd.write(fld)
    log.info(f'Writing the design RosettaScript file: {wpaths["designRS"]}\n')
    with wpaths['designRS'].open('w') as fd:
        fd.write(dsg)

    if TBcore.get_option('system', 'verbose'):
        log.info('Writing weigth 0 patches for folding: {}, {}\n'.format(wpaths['wts0F'], wpaths['wts0_patchF']))
        log.info('Writing weigth 1 patches for folding: {}, {}\n'.format(wpaths['wts1F'], wpaths['wts1_patchF']))
        log.info('Writing weigth 2 patches for folding: {}, {}\n'.format(wpaths['wts2F'], wpaths['wts2_patchF']))
        log.info('Writing weigth 3 patches for folding: {}, {}\n'.format(wpaths['wts3F'], wpaths['wts3_patchF']))
        log.info('Writing weigth 5 patches for folding: {}, {}\n'.format(wpaths['wts5F'], wpaths['wts5_patchF']))
    with wpaths['wts0F'].open('w') as fd: fd.write(wts[0])
    with wpaths['wts0_patchF'].open('w') as fd: fd.write(wts[1])
    with wpaths['wts1F'].open('w') as fd: fd.write(wts[2])
    with wpaths['wts1_patchF'].open('w') as fd: fd.write(wts[3])
    with wpaths['wts2F'].open('w') as fd: fd.write(wts[4])
    with wpaths['wts2_patchF'].open('w') as fd: fd.write(wts[5])
    with wpaths['wts3F'].open('w') as fd: fd.write(wts[6])
    with wpaths['wts3_patchF'].open('w') as fd: fd.write(wts[7])
    with wpaths['wts5F'].open('w') as fd: fd.write(wts[8])
    with wpaths['wts5_patchF'].open('w') as fd: fd.write(wts[9])

    data['script']['folding'] = wpaths['foldRS']
    data['script']['design'] = wpaths['designRS']
    data['cmd']['folding'].append(wpaths['foldRS'])
    data['cmd']['design'].append(wpaths['designRS'])
    return data


def commands( case: Case, nstruct: int, design_nstruct: int, data: Dict, wpaths: Dict ) -> Dict:
    """Create the full commands to execute Rosetta.
    """
    out_prefix = (case.name if not TBcore.get_option('slurm', 'use') else '_'.join([case.name, '${SLURM_ARRAY_TASK_ID}'])) + '_'
    nstruct = nstruct if not TBcore.get_option('slurm', 'use') else math.ceil(nstruct / TBcore.get_option('slurm', 'array'))

    commons = ['-overwrite', '-in:ignore_unrecognized_res', '-in:ignore_waters', '-out:file:silent_struct_type',
               'binary', '-out:mute', 'protocols.abinitio', 'protocols.moves', 'core.optimization']

    flded = str(wpaths['outdir'].joinpath(out_prefix + 'funfol')) + '.silent'
    dsgnd = str(wpaths['outdir'].joinpath(out_prefix + 'des')) + '.silent'
    prefix1 = out_prefix + 'funfol_'
    prefix2 = out_prefix + 'des_'
    data['cmd']['folding'].extend(['-in:file:s', str(wpaths['pdb']), '-out:prefix', prefix1, '-out:file:silent', flded])
    data['cmd']['design'].extend(['-in:file:silent', flded, '-out:prefix', prefix2, '-out:file:silent', dsgnd])
    data['cmd']['folding'].extend(['-nstruct', str(nstruct)])
    data['cmd']['design'].extend(['-nstruct', str(design_nstruct)])
    data['cmd']['folding'].extend(commons)
    data['cmd']['design'].extend(commons)
    return data


def execute( log: Logger, data: Dict, wpaths: Dict ) -> Dict:
    """Run Rosetta.
    """
    if TBcore.get_option('slurm', 'use'):
        slurm_file = wpaths['main'].joinpath('submit_funfoldes.sh')
        log.notice(f'Submission file at {slurm_file}')
        with slurm_file.open('w') as fd:
            fd.write(TButil.slurm_header() + '\n' )
            for k in ['folding', 'design']:
                cmd = ['srun', ]
                cmd.extend(data['cmd'][k])
                fd.write(' '.join([str(x) for x in cmd]) + '\n')

        log.info('Submiting jobs to SLURM... this might take a while\n')
        TButil.submit_slurm(log, slurm_file)
    else:
        for k in ['folding', 'design']:
            log.notice(f'EXECUTE: {" ".join([str(x) for x in data["cmd"][k]])}')
            run([str(x) for x in data['cmd'][k]], stdout=DEVNULL)
    return data


def update_data( log: Logger, data: Dict, wpaths: Dict ) -> Dict:
    """Update data to link final files.
    """
    data['silent_files']['folding'] = list(wpaths['outdir'].glob('*_funfol.silent'))
    data['silent_files']['design'] = list(wpaths['outdir'].glob('*_des.silent'))
    data['minisilent']['folding'] = wpaths['main'].joinpath('output_funfol.minisilent.gz')
    data['minisilent']['design'] = wpaths['main'].joinpath('output_des.minisilent.gz')
    for k in data['minisilent']:
        log.info(f'Generating minisilent file at {data["minisilent"][k]}\n')
        fd = gzip.open( data['minisilent'][k], "wb" )
        for line, _, _, _ in open_rosetta_file([str(x) for x in data['silent_files'][k]], True, check_symmetry=False ):
            fd.write(line.encode('utf-8'))
    return data
