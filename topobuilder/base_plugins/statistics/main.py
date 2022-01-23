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
import glob
import textwrap
from pathlib import Path
from typing import List, Dict, Optional
from subprocess import run
from inspect import getmembers, isfunction

# External Libraries
import pandas as pd

# This Library
from topobuilder.workflow import Node, NodeDataError, NodeMissingError
from topobuilder.case import Case
from .core import core
import topobuilder.core as TBcore
import topobuilder.utils as TButil
from .utils import pdb2fasta, pdb2a3m
from . import orders


__all__ = ['statistics']


class statistics( Node ):
    """Various statistics on the sequence and structure level are computed depending on available scripts.

    .. note::
        Depends on the ``statistic.molprobity`` configuration option.
        Depends on the ``statistic.tmalign`` configuration option.
        Depends on the ``statistic.trrosetta_repo`` configuration option.
        Depends on the ``statistic.trrosetta_wts`` configuration option.
        Depends on the ``statistic.trrosetta_env`` configuration option.

    .. caution::
        In order to execute this :class:`.Node`, we highly recommend to install `trRosetta` with all dependencies.
        The external conda environment can be specified in the ``statistic.trrosetta_env`` configuration option.

    .. admonition:: To Developers

        Due to its use in multiple :class:`.Node`, functions to deal with this :class:`.Node` are mostly located
        in the respective module file and external scripts are locate in this :class:`.Node` directory.

    :param loop_range: Expected loop length is calculated from the euclidian distance between two secondary
        structures. This attribute adds a window of ``loop_range`` residues under and over the calculated
        length.
    :param source: Plugin designs come from, e.g. :class:`funfoldes`.
    :param stage: The type of design, e.g. folded or designed.
    :param analysis: Geometric or quality assessment.
    :param metric: Type of geometric or quality assessment.

    :raises:
        :NodeDataError: On **check**. If the required fields to be executed are not there.
    """
    REQUIRED_FIELDS = ()
    RETURNED_FIELDS = ()
    VERSION = 'v1.0'

    def __init__( self, tag: int,
                        source: str,
                        stage: str,
                        analysis: str,
                        metric: Optional[str] = None,
                        **kwargs ) -> str:
        super(statistics, self).__init__(tag)

        self.source = source
        self.stage = stage
        self.analysis = analysis
        self.metric = metric


    def single_check( self, dummy: Dict ) -> Dict:
        case = Case(dummy)

        # Check what it needs
        for itag in self.REQUIRED_FIELDS:
            if case[itag] is None:
                raise NodeDataError(f'Field "{itag}" is required')

        # Include what keywords it adds (in this instance, nothing)
        if self.analysis == 'geometry':
            case['metadata'].setdefault('statistic', {}).setdefault('geometry', '')
        if self.analysis == 'quality':
            case['metadata'].setdefault('statistic', {}).setdefault('quality', '')
        return case.data

    def single_execute( self, data: Dict ) -> Dict:
        case = Case(data)

        # Generate the folder tree for a single connectivity.
        wfolder = case.connectivities_paths[0].joinpath(f'statistic/{self.source}_{self.stage}/')
        wfolder.mkdir(parents=True, exist_ok=True)
        # Generate internal folder
        thisfolder = wfolder.joinpath('_pdb_files')
        thisfolder.mkdir(parents=True, exist_ok=True)

        # Commands
        commands = []

        # Get data by source
        if len(os.listdir(thisfolder)) == 0:
                commands.extend(self.funfoldes2pdb(case, thisfolder))

        # Load analysis commands
        if self.analysis == 'geometry':
            if not os.path.exists(str(wfolder.joinpath('_geometry.1.csv'))):
                commands.append(self.geometry(case, wfolder, thisfolder))

        if self.analysis == 'quality':
            if not (os.path.exists(str(wfolder.joinpath(f'_{self.metric}.1.csv'))) or
                    os.path.exists(str(wfolder.joinpath(f'_{self.metric}.1.txt'))) ):
                commands.append(self.quality(case, wfolder, thisfolder))

        # Execute
        if commands != []:
            self.execute_commands(commands, wfolder)

        # Postprocess
        self.postprocess(wfolder)

        if self.analysis == 'geometry':
            case['metadata'].setdefault('statistic',
                                        {}).setdefault('geometry',
                                                       str(thisfolder.joinpath('${SLURM_ARRAY_TASK_ID}')))

        if self.analysis == 'quality':
            case['metadata'].setdefault('statistic',
                                        {}).setdefault('quality',
                                                       str(thisfolder.joinpath('${SLURM_ARRAY_TASK_ID}')))

        return case


    def funfoldes2pdb( self, case: Case, wfolder: Path ) -> List:
        """
        """
        if self.stage == 'folding':
            silent_files = case['metadata.funfoldes.silent_files.folding']
        elif self.stage == 'design':
            silent_files = case['metadata.funfoldes.silent_files.design']
        if silent_files is None:
            raise NodeMissingError('There is no output data from the funfoldes plugin.')

        extract_pdb = Path(str(Path(TBcore.get_option('rosetta', 'scripts')).resolve()
                            ).replace('rosetta_scripts.', 'extract_pdbs.'))
        if not extract_pdb.is_file() or not os.access(str(extract_pdb), os.X_OK):
            raise NodeDataError(f'Cannot find executable {extract_pdb}')

        if not TBcore.get_option('slurm', 'use'):
            cmd = [extract_pdb, '-in:file:silent']
            cmd.extend(silent_files)
            cmd.extend(['-out:prefix', str(wfolder) + '/'])
        else:
            indir = str(wfolder.joinpath('${SLURM_ARRAY_TASK_ID}'))
            cmd = ['srun', extract_pdb, '-in:file:silent']
            if self.stage == 'folding':
                cmd.append(os.path.commonprefix([str(x) for x in silent_files]) + '${SLURM_ARRAY_TASK_ID}_funfol.silent')
            elif self.stage == 'design':
                cmd.append(os.path.commonprefix([str(x) for x in silent_files]) + '${SLURM_ARRAY_TASK_ID}_des.silent')
            cmd.extend(['-out:prefix', str(indir) + '/'])
        return [['mkdir', '-p', indir] if TBcore.get_option('slurm', 'use') else '', cmd]


    def hybridize2pdb( self, case: Case, wfolder: Path ) -> List:
        """
        """
        if self.stage == 'assembly':
            silent_files = case['metadata.hybridize.silent_files.assembly']
        elif self.stage == 'design':
            silent_files = case['metadata.hybridize.silent_files.design']
        if silent_files is None:
            raise NodeMissingError('There is no output data from the hybridize plugin.')

        extract_pdb = Path(str(Path(TBcore.get_option('rosetta', 'scripts')).resolve()
                          ).replace('rosetta_scripts.', 'extract_pdbs.'))
        if not extract_pdb.is_file() or not os.access(str(extract_pdb), os.X_OK):
            raise NodeDataError(f'Cannot find executable {extract_pdb}')

        if not TBcore.get_option('slurm', 'use'):
            cmd = [extract_pdb, '-in:file:silent']
            cmd.extend(silent_files)
            cmd.extend(['-out:prefix', str(wfolder) + '/'])
        else:
            indir = str(wfolder.joinpath('${SLURM_ARRAY_TASK_ID}'))
            cmd = ['srun', extract_pdb, '-in:file:silent']
            if self.stage == 'assembly':
                cmd.append(os.path.commonprefix([str(x) for x in silent_files]) + '${SLURM_ARRAY_TASK_ID}_hyb.silent')
            elif self.stage == 'design':
                cmd.append(os.path.commonprefix([str(x) for x in silent_files]) + '${SLURM_ARRAY_TASK_ID}_des.silent')
            cmd.extend(['-out:prefix', str(indir) + '/'])
        return [['mkdir', '-p', indir] if TBcore.get_option('slurm', 'use') else '', cmd]


    def geometry( self, case: Case, wfolder: Path, thisfolder: Path ) -> List:
        """
        """
        cfile = case.write(wfolder.joinpath('current_case'))
        cmd = ['python', str(Path(__file__).parent.joinpath('geometry.py')), '-case', str(cfile), '-indir']
        if not TBcore.get_option('slurm', 'use'):
            cmd.append(str(thisfolder) + '/')
            cmd.extend(['-out', str(wfolder.joinpath('_geometry.1.csv'))])
        else:
            cmd.append(str(thisfolder.joinpath('${SLURM_ARRAY_TASK_ID}')))
            cmd.extend(['-out', str(wfolder.joinpath('_geometry.${SLURM_ARRAY_TASK_ID}.csv'))])
        return cmd


    def quality( self, case: Case, wfolder: Path, thisfolder: Path ) -> List:
        """
        """
        if self.metric == 'molprobity':
            self.log.notice(f'Quality assessment of decoys: {self.metric}')
            cfile = case.write(wfolder.joinpath('current_case'))
            cmd = [f'{core.get_option("statistics", "molprobity")}']
            if not TBcore.get_option('slurm', 'use'):
                cmd.append(str(thisfolder) + '/')
                cmd.extend(['> ' + str(wfolder.joinpath('_molprobity.1.txt'))])
            else:
                cmd.append(str(thisfolder.joinpath('${SLURM_ARRAY_TASK_ID}')))
                cmd.extend(['> ' + str(wfolder.joinpath('_molprobity.${SLURM_ARRAY_TASK_ID}.txt'))])

        if self.metric == 'proq4':
            self.log.notice(f'Quality assessment of decoys: {self.metric}')
            cmd = [f'python {Path(__file__).parent.joinpath("calc_proq4.py")}']
            if not TBcore.get_option('slurm', 'use'):
                cmd.extend(['-p', '1', '-f', str(thisfolder), '-o', str(wfolder)])
            else:
                cmd.extend(['-p', '${SLURM_ARRAY_TASK_ID}', '-f', str(thisfolder.joinpath('${SLURM_ARRAY_TASK_ID}')), '-o', str(wfolder)])

        if self.metric == 'trRosetta':
            self.log.notice(f'Quality assessment of decoys: {self.metric}')
            cfile = case.write(wfolder.joinpath('current_case'))
            cmd = [f'source {core.get_option("statistics", "trrosetta_env")}']
            if not TBcore.get_option('slurm', 'use'):
                # Predict restraints and build model for each
                for pdbfile in glob.iglob(str(thisfolder) + '/*/*.pdb'):
                    a3mfile = pdbfile.replace('.pdb', '.a3m')
                    pdb2a3m(pdbfile, a3mfile)
                    fastafile = pdbfile.replace('.pdb', '.fasta')
                    pdb2fasta(pdbfile, fastafile)
                    base = a3mfile.replace('.a3m', '')
                    slurmarray = a3mfile.split('/')[-2]
                    modelnumber = a3mfile.split('_')[-1].replace('.a3m', '')
                    cmd.extend([f'\npython {core.get_option("statistics", "trrosetta_repo")}network/predict_many.py -m {core.get_option("statistics", "trrosetta_wts")}',
                                str(thisfolder) + '/' + ' ' + str(thisfolder) + '/'])
                    cmd.extend(['\npython', f'{core.get_option("statistics", "trrosetta_repo")}trRosetta_modelling/trRosetta.py {base}.npz {base}.a3m {base}_tr001.1.pdb'])
                    if TBcore.get_option('tmalign', 'script'):
                        cmd.extend([f'\n{core.get_option("statistics", "tmalign")} {base}.pdb {base}_tr001.1.pdb', '>',
                                    f'{str(thisfolder)}/_trRosetta.{modelnumber}.1'])
                        cmd.extend(['\npython', Path(__file__).parent.joinpath('parse_tmalign.py'), '-f', f'{str(thisfolder)}/',
                                    '-o', f'{str(thisfolder)}/_trRosetta.1.csv'])
            else:
                # Predict restraints
                with open(f'{wfolder}/exec_trRosetta_predict.sh', 'w') as f:
                    f.write(f'source {core.get_option("statistics", "trrosetta_env")}\n\n')
                    f.write('SLURM_ARRAY_TASK_ID=$1\n\n')
                    slurmarrays = []
                    for pdbfile in glob.iglob(str(thisfolder) + '/*/*.pdb'):
                        a3mfile = pdbfile.replace('.pdb', '.a3m')
                        pdb2a3m(pdbfile, a3mfile)
                        fastafile = pdbfile.replace('.pdb', '.fasta')
                        pdb2fasta(pdbfile, fastafile)
                        base = a3mfile.replace('.a3m', '')
                        slurmarray = a3mfile.split('/')[-2]
                        modelnumber = a3mfile.split('_')[-1].replace('.a3m', '')
                        slurmarrays.append(slurmarray)
                    for slurmarray in list(set(slurmarrays)):
                        f.write(f'\n\nif (( ${{SLURM_ARRAY_TASK_ID}} == {slurmarray} )); then')
                        f.write(f'\npython {core.get_option("statistics", "trrosetta_repo")}network/predict_many.py -m {core.get_option("statistics", "trrosetta_wts")}' + ' ' + str(thisfolder.joinpath(slurmarray)) + ' ' + str(thisfolder.joinpath(slurmarray)))
                        f.write('\nfi')
                cmd.extend([f'\nbash {wfolder}/exec_trRosetta_predict.sh ${{SLURM_ARRAY_TASK_ID}}'])

                # Build model for each
                with open(f'{wfolder}/exec_trRosetta_relax.sh', 'w') as f:
                    f.write(f'source {core.get_option("statistics", "trrosetta_env")}\n\n')
                    f.write('SLURM_ARRAY_TASK_ID=$1\n\n')
                    for pdbfile in glob.iglob(str(thisfolder) + '/*/*.pdb'):
                        base = pdbfile.replace('.pdb', '')
                        slurmarray = pdbfile.split('/')[-2]
                        modelnumber = pdbfile.split('_')[-1].replace('.pdb', '')
                        f.write(f'\n\nif (( ${{SLURM_ARRAY_TASK_ID}} == {slurmarray} )); then')
                        f.write(f'\npython {core.get_option("statistics", "trrosetta_repo")}trRosetta_modelling/trRosetta.py {base}.npz {base}.a3m {base}_tr001.{slurmarray}.pdb')
                        f.write(f'\n{core.get_option("statistics", "tmalign")} {base}.pdb {base}_tr001.{slurmarray}.pdb > {str(thisfolder.joinpath(slurmarray))}/_trRosetta.{modelnumber}.{slurmarray}')
                        f.write('\nfi')
                cmd.extend([f'\nbash {wfolder}/exec_trRosetta_relax.sh ${{SLURM_ARRAY_TASK_ID}}'])

                # Parse TM aligns
                with open(f'{wfolder}/exec_trRosetta_TMalign.sh', 'w') as f:
                    f.write(f'source {core.get_option("statistics", "trrosetta_env")}\n\n')
                    for slurmarray in list(set(slurmarrays)):
                        if core.get_option('statistics', 'tmalign'):
                            f.write(f'\n\nif (( ${{SLURM_ARRAY_TASK_ID}} == {slurmarray} )); then')
                            f.write('\npython ' + str(Path(__file__).parent.joinpath('parse_tmalign.py')) + ' -f ' + f'{str(thisfolder.joinpath(slurmarray))}/' + ' -o ' + f'{str(wfolder)}/_trRosetta.{slurmarray}.csv')
                            f.write('\nfi')
                cmd.extend([f'\nbash {wfolder}/exec_trRosetta_TMalign.sh ${{SLURM_ARRAY_TASK_ID}}'])

                # # Predict restraints and build model for each
                # for pdbfile in glob.iglob(str(thisfolder) + '/*/*.pdb'):
                #     a3mfile = pdbfile.replace('.pdb', '.a3m')
                #     pdb2a3m(pdbfile, a3mfile)
                #     fastafile = pdbfile.replace('.pdb', '.fasta')
                #     pdb2fasta(pdbfile, fastafile)
                #     base = a3mfile.replace('.a3m', '')
                #     slurmarray = a3{wfolder}/mfile.split('/')[-2]
                #     modelnumber = a3mfile.split('_')[-1].replace('.a3m', '')
                #     cmd.extend([f'\n\nif (( ${{SLURM_ARRAY_TASK_ID}} == {slurmarray} )); then'])
                #     cmd.extend([f'\npython {core.get_option("statistics", "trrosetta_repo")}network/predict_many.py -m {core.get_option("statistics", "trrosetta_wts")}',
                #                 str(thisfolder.joinpath(slurmarray)) + ' ' + str(thisfolder.joinpath(slurmarray))])
                #     cmd.extend(['\npython', f'{core.get_option("statistics", "trrosetta_repo")}trRosetta_modelling/trRosetta.py {base}.npz {base}.a3m {base}_tr001.{slurmarray}.pdb'])
                #     if core.get_option('statistics', 'tmalign'):
                #         cmd.extend([f'\n{core.get_option("statistics", "tmalign")} {base}.pdb {base}_tr001.{slurmarray}.pdb', '>',
                #                     f'{str(thisfolder.joinpath(slurmarray))}/_trRosetta.{modelnumber}.{slurmarray}'])
                #         cmd.extend(['\npython', Path(__file__).parent.joinpath('parse_tmalign.py'), '-f', f'{str(thisfolder.joinpath(slurmarray))}/',
                #                     '-o', f'{str(wfolder)}/_trRosetta.{slurmarray}.csv'])
                #     cmd.extend(['\nfi'])
        return cmd


    def execute_commands( self, cmd: List, wfolder: Path ):
        """
        """
        if not TBcore.get_option('slurm', 'use'):
            for c in list(filter(None, cmd)):
                run(c)
        else:
            if self.metric:
                slurm_file = wfolder.joinpath(f'submit_analytics_{self.analysis}_{self.metric}.sh')
            else:
                slurm_file = wfolder.joinpath(f'submit_analytics_{self.analysis}.sh')
            with slurm_file.open('w') as fd:
                fd.write(TButil.slurm_header() + '\n')
                fd.write(TButil.slurm_pyenv() + '\n')
                for c in cmd:
                    fd.write(' '.join([str(x) for x in c]) + '\n')
            TButil.submit_slurm(self.log, slurm_file)


    def postprocess( self, wfolder: Path ):
        """
        """
        if self.analysis == 'geometry':
            df = pd.concat([pd.read_csv(x) for x in Path(wfolder).glob('_geometry.*.csv')])
            df.to_csv(wfolder.joinpath('geometry.csv'), index=False)
        if self.analysis == 'quality':
            # MolProbity
            if os.path.exists(wfolder.joinpath('_molprobity.1.txt')):
                _getter_ = ['description', 'MolProbityScore', 'Mol_pct_rank', 'clashscore',
                            'pct_rank', 'cbeta>0.25', 'ramaOutlier', 'ramaAllowed', 'ramaFavored']
                df = pd.concat([pd.read_csv(x, sep=":", skiprows=2).rename(columns={'#pdbFileName': 'description'})[_getter_]
                                    for x in Path(wfolder).glob('_molprobity.*.txt')])
                df.to_csv(wfolder.joinpath('molprobity.csv'), index=False)

            # ProQ4
            if os.path.exists(wfolder.joinpath('_proq4.1.csv')):
                df = pd.concat([pd.read_csv(x) for x in Path(wfolder).glob('_proq4.*.csv')])
                df.to_csv(wfolder.joinpath('proq4.csv'), index=False)

            # trRosetta
            if os.path.exists(wfolder.joinpath('_trRosetta.1.csv')):
                df = pd.concat([pd.read_csv(x) for x in Path(wfolder).glob('_trRosetta.*.csv')])
                df.to_csv(wfolder.joinpath('trRosetta.csv'), index=False)
