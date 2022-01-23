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
from logbook import Logger
from topobuilder.case import Case
import topobuilder.core as TBcore
import topobuilder.utils as TButil


def make_folder_structure( case: Case, statmode: str ) -> Dict:
    """
    """

    # Generate the folder tree for a single connectivity.
    wfolder = case.connectivities_paths[0].joinpath('statistic')
    wfolder.mkdir(parents=True, exist_ok=True)

    # Generate internal folder
    thisfolder = wfolder.joinpath(statmode)
    thisfolder.mkdir(parents=True, exist_ok=True)

    return {'global': wfolder, 'main': thisfolder}


def make_directed_sketch( log: Logger, case: Case, folder: Path ) -> Path:
    """
    """

    pdbfile = folder.joinpath('directed_sketch.pdb')
    structure, _ = TButil.build_pdb_object(case.apply_topologies()[0].ordered_structures, 3)
    log.notice(f'New File: Writing structure {pdbfile}')
    structure.write(output_file=str(pdbfile), format='pdb', clean=True, force=True)
    return pdbfile


def count_single_master_matches( log: Logger, pdbfile: Path, folder: Path ) -> pd.DataFrame:
    """
    """

    createpds = TButil.createPDS(pdbfile)
    log.notice(f'EXECUTE: {createpds}')
    run(createpds, stdout=DEVNULL)
    masters = TButil.master_best_each(pdbfile.with_suffix('.pds'), folder.joinpath('_master'), 5)

    unimaster = folder.joinpath('match.master')

    if not TBcore.get_option('slurm', 'use'):
        no_slurm(cmd, current_case_file, current_sse, unimaster, imaster, unidata.with_suffix(''))
    else:
        with_slurm(cmd, current_case_file, current_sse, unimaster, imaster, unidata)

    return {'matches': unimaster, 'stats': unidata, 'corrections': None,
            'layers': list(set([x[0] for x in current_sse.split('.')]))}


    data = submit_searches(masters, stepfolder, current_case_file, '.'.join([x['id'] for x in sses]))


def pdb2a3m(infile, outfile):
    """
    Creates a m3a file from a pdb file.
    """
    letters = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H',
               'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
               'TYR':'Y','VAL':'V'}

    with open(infile, "r") as infl:
        with open(outfile, "w") as otfl:
            prev = '-1'
            for line in infl:
                toks = line.split()
                if len(toks)<1: continue
                if toks[0] != 'ATOM': continue
                if toks[2] == 'CA':
                    otfl.write('%c' % letters[toks[3]])

            otfl.write('\n')


def pdb2fasta(infile, outfile):
    """
    Creates a fasta file from a pdb file.
    """
    letters = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H',
               'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
               'TYR':'Y','VAL':'V'}

    with open(infile, "r") as infl:
        base = os.path.basename(infile).replace('.pdb', '')
        with open(outfile, "w") as otfl:
            prev = '-1'
            otfl.write(f'> {base}\n')
            for line in infl:
                toks = line.split()
                if len(toks)<1: continue
                if toks[0] != 'ATOM': continue
                if toks[2] == 'CA':
                    otfl.write('%c' % letters[toks[3]])

            otfl.write('\n')
