# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>


"""
# Standard Libraries
import sys
import os
import shlex
from pathlib import Path
from typing import Optional, Tuple, List, Union
from ast import literal_eval
from tempfile import NamedTemporaryFile

# External Libraries
import pandas as pd
import numpy as np
from logbook import Logger

# This Library
import topobuilder.core as TBcore
from topobuilder.workflow import NodeDataError


def get_master_exes() -> Tuple[str, str]:
    """Provide the path to the ``MASTER`` executables.

    .. note::
        Depends on the ``master.master`` configuration option
        Depends on the ``master.create`` configuration option

    :return: The ``master`` and ``createPDS`` executable paths. In this order.
    """
    exes = []
    for name, option_name in [('master', 'master'), ('createPDS', 'create')]:
        # Get the executable and check that it works
        exe = TBcore.get_option('master', option_name)
        if not os.access(exe, os.X_OK):
            raise IOError('Unable to find a proper executable for {} at {}'.format(name, exe))
        exes.append(exe)
    return exes


def pds_database( log: Logger,
                  filter: Optional[Union[str, Path]] = None,
                  ) -> Tuple[Path, List]:
    """Provide the list of target PDS as a file and a list.

    .. note::
        Depends on the ``master.pds`` configuration option.

    :param log: Job logger.
    :param filter: File containting the target subset, otherwise all PDS database is considered.
    """
    # @TODO: PDS-FILTER
    pds_file = Path(TBcore.get_option('master', 'pds'))
    if pds_file.is_file():
        pds_list = [line.strip() for line in open(pds_file).readlines() if len(line.strip()) > 0]
    elif pds_file.is_dir():
        pds_list = [str(x.resolve()) for x in pds_file.glob('*/*.pds')]
    else:
        raise NodeDataError('The provided MASTER database directory/list file cannot be found.')

    # Even if a PDS file already exists, we always create a temporary file so that we can
    # manage different versions of the PDS database in different Nodes.
    f = NamedTemporaryFile(mode='w', delete=False)
    log.info(f'Temporary file for PDS database: {f.name}')
    [f.write(x + '\n') for x in pds_list]
    f.close()
    pds_file = Path(f.name)
    return pds_file, pds_list


def createPDS( infile: Union[Path, str], outfile: Optional[str] = None ) -> List[str]:
    """Make the createPDS command call.

    .. note::
        Depends on the ``master.create`` configuration option.

    :param infile: PDB file to convert.
    :param outfile: Name of the expected PDS output.
    """
    _, createPDS = get_master_exes()
    infile = Path(infile)
    if not infile.is_file():
        raise NodeDataError(f'Unable to find structure file {infile}')
    outfile = outfile if outfile is not None else infile.with_suffix('.pds')
    return shlex.split(f'{createPDS} --type query --pdb {str(infile)} --pds {str(outfile)}  --dCut 100')


def master_best_each( log: Logger,
                      infile: Union[Path, str],
                      outdir: Union[Path, str],
                      rmsd: Optional[float] = 5.0
                      ) -> List[List[str]]:
    """Create one MASTER call for each PDS file with the --topN 1 flag.

    .. note::
        Depends on the ``master.master`` configuration option
        Depends on the ``master.pds`` configuration option
    """
    master, _ = get_master_exes()
    _, pds_list = pds_database(log)

    infile = Path(infile)
    if not infile.is_file():
        raise IOError('Unable to find PDS file {}'.format(infile))
    outdir = Path(outdir)
    if not outdir.is_dir():
        outdir.mkdir(parents=True, exist_ok=True)
    createbash = '{0} --query {1} --target {2} --rmsdCut {3}  --topN 1 --matchOut {4}'

    cmds = []
    for pds in pds_list:
        tid = str(Path(Path(pds).name).with_suffix(''))
        outfile = outdir.joinpath('{}.master'.format(tid))
        cmds.append(shlex.split(createbash.format(master, infile, pds, rmsd, outfile)))
    return cmds


def master_fixedgap( query: Path,
                     pds_list: Path,
                     master_out: Path,
                     mdis: int,
                     Mdis: int,
                     rmsd_cut: float
                     ) -> List[List[str]]:
    """Create the :term:`MASTER` executable for a fixed residue distance between matches.

    :param query: Query PDS file.
    :param pds_list: File with the database list of PDS targets.
    :param master_out: Output file with the MASTER results.
    :param mdis: Minimum distance.
    :param Mdis: Maximum distance.
    :param rmsd_cut: RMSD search cutoff.
    """
    master, _ = get_master_exes()
    masterbash = f'{master} --query {str(query)} --targetList {str(pds_list)} --rmsdCut {rmsd_cut} '
    masterbash += f'--matchOut {str(master_out)} --gapLen {mdis}-{Mdis}'
    return shlex.split(masterbash)


def master_groupedgap( query: Path,
                       pds_list: Path,
                       master_out: Path,
                       gaps: str,
                       rmsd_cut: float
                       ) -> List[List[str]]:
    """Create the :term:`MASTER` executable for a fixed residue distance between matches.

    :param query: Query PDS file.
    :param pds_list: File with the database list of PDS targets.
    :param master_out: Output file with the MASTER results.
    :param gaps: Allowed lengths between the gaps.
    :param rmsd_cut: RMSD search cutoff.
    """
    master, _ = get_master_exes()
    masterbash = f'{master} --query {str(query)} --targetList {str(pds_list)} --rmsdCut {rmsd_cut} '
    masterbash += f'--matchOut {str(master_out)} --gapLen {gaps}'
    return shlex.split(masterbash)


def master_nogap( query: Path,
                  pds_list: Path,
                  master_out: Path,
                  rmsd_cut: float
                  ) -> List[List[str]]:
    """Create the :term:`MASTER` executable for a fixed residue distance between matches.

    :param query: Query PDS file.
    :param pds_list: File with the database list of PDS targets.
    :param master_out: Output file with the MASTER results.
    :param mdis: Minimum distance.
    :param Mdis: Maximum distance.
    :param rmsd_cut: RMSD search cutoff.
    """
    master, _ = get_master_exes()
    masterbash = f'{master} --query {str(query)} --targetList {str(pds_list)} --rmsdCut {rmsd_cut} '
    masterbash += f'--matchOut {str(master_out)}'
    return shlex.split(masterbash)
