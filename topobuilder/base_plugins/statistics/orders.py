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
import sys
from pathlib import Path
from typing import List, Dict
from subprocess import run, DEVNULL

# External Libraries
import pandas as pd

# This Library
from topobuilder.case import Case
import topobuilder.core as TBcore
import topobuilder.utils as TButil
from . import utils


def metadata( order: str ) -> Dict:
    """Plugin description.

    It includes:

    - ``name``: The plugin identifier.
    - ``Itags``: The metadata tags neccessary to execute.
    - ``Otags``: The metadata tags generated after a successful execution.
    - ``Isngl``: When :data:`True`, input requires single connectivity.
    - ``Osngl``: When :data:`True`, output guarantees single connectivity.
    """
    if order == 'sketch_2_master':
        return {'name': 'statistics:{}'.format(order),
                'Itags': [],
                'Otags': ['sketch_2_master'],
                'Isngl': True,
                'Osngl': True}


@TButil.plugin_conditions(metadata('sketch_2_master'))
def sketch_2_master( case: Case, **kwargs ):
    """
    """
    # params
    rmsd = kwargs.get('rmsd', 5)

    # Get working paths
    wpath = utils.make_folder_structure(case, 'sketch_2_master')

    # Generate structure query
    pdbfile = utils.make_directed_sketch(case, wpath['main'])

    # MASTER search
    createpds = TButil.createPDS(pdbfile)
    #TButil.plugin_bash(createpds)
    run(createpds, stdout=DEVNULL)
    masters = TButil.master_best_each(pdbfile.with_suffix('.pds'), wpath['main'].joinpath('_master'), rmsd)
    data = submit_searches(masters, stepfolder, current_case_file, '.'.join([x['id'] for x in sses]))


    kase.data['metadata']['imaster'].setdefault('step{:02d}'.format(i + 1), data)
    TButil.checkpoint_out(checkpoint, data)
    kase.data['metadata']['corrections'].append(data['corrections'])
    done_l.update(data['layers'])
    corrections = data['corrections']

#
# def with_slurm( cmd: List[str],
#                 current_case_file: Path,
#                 current_sse: str,
#                 unimaster: Path,
#                 imaster: Path,
#                 unidata: Path ):
#     """
#     """
#     # Make bashfile
#     bashcont = []
#     createbash = 'python {0} -case {1} -master {2} -present {3} -out {4}'
#     parts = math.ceil(len(cmd) / TBcore.get_option('slurm', 'array'))
#
#     wwd = unimaster.parent.parent
#     cwd = Path().cwd()
#     os.chdir(str(wwd))
#
#     for i, com in enumerate(cmd):
#         cmd[i][2] = str(Path(com[2]).relative_to(wwd))
#         cmd[i][-1] = str(Path(com[-1]).relative_to(wwd))
#
#     for ii, cp in enumerate(cmd):
#         cmd[ii][-1] = cp[-1] + '_${SLURM_ARRAY_TASK_ID}'
#     for j, i in enumerate(range(0, len(cmd), parts)):
#         sumfile = unimaster.parent.joinpath('_${SLURM_ARRAY_TASK_ID}.master').relative_to(wwd)
#         datfile = unimaster.parent.joinpath('_${SLURM_ARRAY_TASK_ID}.geo').relative_to(wwd)
#         bashcont.append('if (( ${{SLURM_ARRAY_TASK_ID}} == {} )); then'.format(j + 1))
#         bashcont.extend([' '.join(x) for x in cmd[i:i + parts]])
#         bashcont.append('cat {0} > {1}'.format(Path(cmd[-1][-1]).parent.joinpath('*_${SLURM_ARRAY_TASK_ID}'), sumfile))
#         bashcont.append(createbash.format(imaster, current_case_file.relative_to(wwd),
#                                           sumfile, current_sse, datfile))
#         bashcont.append('fi')
#     with unimaster.parent.joinpath('submit.sh').relative_to(wwd).open('w') as fd:
#         fd.write(TButil.slurm_header())
#         fd.write(TButil.slurm_pyenv())
#         fd.write('\n'.join(bashcont))
#
#     TButil.submit_slurm(unimaster.parent.joinpath('submit.sh').relative_to(wwd))
#     TButil.plugin_filemaker('Creating geometric coordinate file {}'.format(unidata))
#     allCSV = [str(x) for x in unimaster.parent.relative_to(wwd).glob('_*.geo.csv')]
#     pd.concat([pd.read_csv(x) for x in allCSV]).to_csv(unidata.relative_to(wwd), index=False)
#     TButil.plugin_filemaker('Creating MASTER search file {}'.format(unimaster))
#     with unimaster.relative_to(wwd).open('w') as fd:
#         for x in unimaster.parent.glob('_*.master'):
#             with x.relative_to(wwd).open() as fi:
#                 fd.write(''.join(fi.readlines()))
#     os.chdir(str(cwd))
#
#
# def no_slurm( cmd: List[str],
#               current_case_file: Path,
#               current_sse: str,
#               unimaster: Path,
#               imaster: Path,
#               unidata: Path ):
#     """
#     """
#     # Search on MASTER
#     result = []
#     for com in cmd:
#         TButil.plugin_bash(com)
#         run(com, stdout=DEVNULL)
#         outf = Path(com[-1])
#         if outf.is_file():
#             result.append(str(outf))
#     result.insert(0, 'cat')
#     TButil.plugin_filemaker('Unify matches at {0}'.format(unimaster))
#     with unimaster.open('w') as fd:
#         run(result, stdout=fd)
#     result[0] = 'rm'
#     run(result[:-2], stdout=DEVNULL)
#
#     # Analyze
#     createbash = 'python {0} -case {1} -master {2} -present {3} -out {4}'
#     cmd = shlex.split(createbash.format(imaster, current_case_file, unimaster, current_sse, unidata))
#     TButil.plugin_bash(cmd)
#     run(cmd, stdout=DEVNULL)
