# -*- coding: utf-8 -*-
"""
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
from typing import Union, Dict, Tuple
import sys
import textwrap
from subprocess import run, DEVNULL

# External Libraries
from rstoolbox.io import open_rosetta_file

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
    folders = case.connectivities_paths[0].joinpath('hybridize')
    folders.mkdir(parents=True, exist_ok=True)
    outdir = folders.joinpath('outputs')
    outdir.mkdir(parents=True, exist_ok=True)
    pdb_file = folders.joinpath('template_sketch.pdb')
    hyz_assbly_file = folders.joinpath('hybridize_assembly.xml')
    hyz_design_file = folders.joinpath('hybridize_design.xml')
    checkpoint = folders.joinpath('checkpoint.json')

    return {'main': folders,                # Main plugin folder
            'outdir': outdir,               # Folder to produce the Rosetta outputs
            'pdb': pdb_file,                # Path to the template PDB file
            'assemblyRS': hyz_assbly_file,    # Path to the Assembly script
            'designRS': hyz_design_file,    # Path to the Design script
            'checkpoint': checkpoint        # Path to the Checkpoint file (.json)
            }


def build_template_sketch( log: Logger, case: Case, pdb_file: Union[Path, str] ):
    """Generate the PDB file used as template to fold.
    """
    # 1.1 Make sure the loop_length are specified.
    loop_lengths = case['metadata.loop_lengths']
    looplist = ', '.join([str(x) for x in loop_lengths])
    log.info(f'Gathering a total of {len(loop_lengths)} loops of length: {looplist}\n')
    log.info(f'To apply over {len(case)} secondary structures: {case.connectivities_str[0]}\n')

    # Get the structure.
    node = getattr(TBplugins.source.load_plugin('builder'), 'builder', None)(connectivity=True, pick_aa='V', tag=0)
    case = Case(node.single_execute(case.data))
    sse_list = case.ordered_structures

    pdb, _ = build_pdb_object(log, sse_list, loop_lengths)
    pdb.write(str(pdb_file), format='pdb', clean=True, force=TBcore.get_option('system', 'overwrite'))


def make_scripts( log: Logger,
                  case: Case,
                  wpaths: Dict,
                  data: Dict,
                  natbias: float = 2.5,
                  layer_design: bool = True
                  ) -> Tuple[str, str]:
    """Create the assembly and design scripts.
    """
    fld = TButil.rosettascript(TButil.hybridize(case, wpaths['pdb'], natbias))
    dsg = TButil.rosettascript(TButil.constraint_design(case, natbias, layer_design))

    if TBcore.get_option('system', 'jupyter'):
        ifold = os.getenv('TB_HYBRIDIZE_ASMB_FILE', None)
        idsgn = os.getenv('TB_HYBRIDIZE_DSGN_FILE', None)

        if ifold is None:
            print('-' * 80)
            print('\n\nExpected ASSEMBLY script will be:\n')
            print(fld)
            print(textwrap.dedent("""\n\n If a different script wants to be provided, do so by
    assigning a RosettaScript to os.environ['TB_HYBRIDIZE_ASMB_FILE'] or assign '' to
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
    assigning a RosettaScript to os.environ['TB_HYBRIDIZE_DSGN_FILE'] or assign '' to
    it if you are ok with the default script."""))
        elif idsgn:
            idsgn = Path(idsgn)
            if not idsgn.is_file():
                raise IOError(f'Unknown file {idsgn}')
            dsg = ''.join(list(idsgn.open().readlines()))

        if ifold is None or idsgn is None:
            TButil.exit()

    log.info(f'Writing the assembly RosettaScript file: {wpaths["assemblyRS"]}\n')
    with wpaths['assemblyRS'].open('w') as fd:
        fd.write(fld)
    log.info(f'Writing the design RosettaScript file: {wpaths["designRS"]}\n')
    with wpaths['designRS'].open('w') as fd:
        fd.write(dsg)

    data['script']['assembly'] = wpaths['assemblyRS']
    data['script']['design'] = wpaths['designRS']
    data['cmd']['assembly'].append(wpaths['assemblyRS'])
    data['cmd']['design'].append(wpaths['designRS'])
    return data


def commands( case: Case, nstruct: int, data: Dict, wpaths: Dict ) -> Dict:
    """Create the full commands to execute Rosetta.
    """
    out_prefix = (case.name if not TBcore.get_option('slurm', 'use') else '_'.join([case.name, '${SLURM_ARRAY_TASK_ID}'])) + '_'
    nstruct = nstruct if not TBcore.get_option('slurm', 'use') else math.ceil(nstruct / TBcore.get_option('slurm', 'array'))

    commons = ['-overwrite', '-in:ignore_unrecognized_res', '-in:ignore_waters', '-out:file:silent_struct_type',
               'binary', '-out:mute', 'protocols.abinitio', 'protocols.moves', 'core.optimization']

    flded = str(wpaths['outdir'].joinpath(out_prefix + 'hyb')) + '.silent'
    dsgnd = str(wpaths['outdir'].joinpath(out_prefix + 'des')) + '.silent'
    prefix1 = out_prefix + 'hyb_'
    prefix2 = out_prefix + 'des_'
    data['cmd']['assembly'].extend(['-in:file:s', str(wpaths['pdb']), '-out:prefix', prefix1, '-out:file:silent', flded])
    data['cmd']['design'].extend(['-in:file:silent', flded, '-out:prefix', prefix2, '-out:file:silent', dsgnd])
    data['cmd']['assembly'].extend(['-nstruct', str(nstruct)])
    data['cmd']['design'].extend(['-nstruct', str(10)])
    data['cmd']['assembly'].extend(commons)
    data['cmd']['design'].extend(commons)
    return data


def execute( log: Logger, data: Dict, wpaths: Dict ) -> Dict:
    """Run Rosetta.
    """
    if TBcore.get_option('slurm', 'use'):
        slurm_file = wpaths['main'].joinpath('submit_hybridize.sh')
        TButil.plugin_filemaker('Submission file at {}'.format(slurm_file))
        with slurm_file.open('w') as fd:
            fd.write(TButil.slurm_header() + '\n' )
            for k in ['assembly', 'design']:
                cmd = ['srun', ]
                cmd.extend(data['cmd'][k])
                fd.write(' '.join([str(x) for x in cmd]) + '\n')

        if TBcore.get_option('system', 'verbose'):
            sys.stdout.write('Submiting jobs to SLURM... this might take a while\n')
        TButil.submit_slurm(slurm_file)
    else:
        for k in ['assembly', 'design']:
            log.notice(f'EXECTUE: {" ".join([str(x) for x in data["cmd"][k]])}')
            run([str(x) for x in data['cmd'][k]], stdout=DEVNULL)
    return data


def update_data( log: Logger, data: Dict, wpaths: Dict ) -> Dict:
    """Update data to link final files.
    """
    data['silent_files']['assembly'] = list(wpaths['outdir'].glob('*_hyb.silent'))
    data['silent_files']['design'] = list(wpaths['outdir'].glob('*_des.silent'))
    data['minisilent']['assembly'] = wpaths['main'].joinpath('output_hyb.minisilent.gz')
    data['minisilent']['design'] = wpaths['main'].joinpath('output_des.minisilent.gz')
    for k in data['minisilent']:
        log.info(f'Generating minisilent file at {data["minisilent"][k]}\n')
        fd = gzip.open( data['minisilent'][k], "wb" )
        for line, _, _, _ in open_rosetta_file([str(x) for x in data['silent_files'][k]], True, check_symmetry=False ):
            fd.write(line.encode('utf-8'))
    return data
