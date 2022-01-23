# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>
.. codeauthor:: Zander Harteveld <zandermilanh@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import Optional, List, Tuple
import os
import sys

# External Libraries
import SBI.databases as SBIdb
import SBI.core as SBIcr
import SBI.structure as SBIstr
import pandas as pd
import numpy as np

# This Library
import topobuilder.core as TBcore
import topobuilder.utils as TButil


__all__ = ['get_steps', 'process_master_geometries']


def get_steps( blist: List[bool] ) -> List[Tuple[int]]:
    """
    """
    def cstp(seq):
        return [seq[i] for i in range(len(seq)) if len(seq[i]) == 1 or (seq[i][0] + 1 == seq[i][1])]

    def f7(seq):
        seen = set()
        seen_add = seen.add
        return [x for x in seq if not (x in seen or seen_add(x))]
    steps = []
    plist = list(range(len(blist)))
    betas = list(np.where(np.asarray(blist))[0])
    if len(betas) == 0:
        steps.append((0, ))
        for i in range(0, len(plist) - 1):
            steps.append(tuple(plist[i: i + 2]))
    else:
        steps.append((betas[0], ))
        for i in range(1, len(betas)):
            if betas[i] > betas[i - 1] + 1:
                steps.append((betas[i], ))
        for i in range(0, len(betas) - 1):
            steps.append(tuple(betas[i: i + 2]))
        for i in range(0, len(plist) - 1):
            steps.append(tuple(plist[i: i + 2]))
        steps = f7(cstp(steps))
    return steps


# assert get_steps([True, False]) == [(0,), (0, 1)]
# assert get_steps([False, True]) == [(1,), (0, 1)]
# assert get_steps([False, False]) == [(0,), (0, 1)]
# assert get_steps([False, False, False]) == [(0,), (0, 1), (1, 2)]
# assert get_steps([False, True, False]) == [(1,), (0, 1), (1, 2)]
# assert get_steps([False, True, True, False]) == [(1,), (1, 2), (0, 1), (2, 3)]
# assert get_steps([False, True, False, True]) == [(1,), (3,), (0, 1), (1, 2), (2, 3)]
# assert get_steps([False, True, True, False, True]) == [(1,), (4,), (1, 2), (0, 1), (2, 3), (3, 4)]


def process_master_geometries( masterdf: pd.DataFrame,
                               structures: List[str],
                               flip: List[str]
                               ) -> pd.DataFrame:
    """Use the PDS matches provided by MASTER to obtain the geometric properties.

    Column modifications:

    =========== ========================================= =======
    column name description                               status
    =========== ========================================= =======
    pds_path    Path to the PDS file containing the match dropped
    pdb_path    Path to the PDB file containing the match new!
    sse         Identifier of the query secondary         new!
                structre.
    layer       Layer against which the data SSE is       new!
                evaluated
    vectors     List of fitting vectors to the SSE.       new!
    planes      List of :class:`tuple` with fitting       new!
                planes and edges to the SSE.
    angles      List of angles between vectors and planes new!
    points      List of distances between the center of   new!
                each vector and the planes
    =========== ========================================= =======
    """
    # Define PDB database.
    with SBIcr.on_option_value('structure', 'format', 'pdb'):
        pdbdb = SBIdb.PDBLink(TBcore.get_option('master', 'pdb'))
        if TBcore.get_option('system', 'verbose'):
            sys.stdout.write('Set PDB database as: {}\n'.format(TBcore.get_option('master', 'pdb')))

    # Prepare data: Sort by structure-chain to avoid multi-reading
    newdf = masterdf.sort_values(['pdb', 'chain'])
    pdb3d = None

    def execute( row, structures, flip ):
        # 1. Download file
        nonlocal pdbdb
        nonlocal pdb3d
        filename, pdb3d = download_pdb(pdbdb, pdb3d, row['pdb'], row['chain'])
        pdb3d = pdb3d['AtomType:CA']
        df = TButil.pdb_geometry_from_rules(pdb3d, list(zip(structures, row['match'], flip)))
        df = df.assign(pdb=[row['pdb'], ] * df.shape[0])
        df = df.assign(chain=[row['chain'], ] * df.shape[0])
        df = df.assign(rmsd=[row['rmsd'], ] * df.shape[0])
        df = df.assign(match=[row['match'], ] * df.shape[0])
        df['pdb_path'] = [filename, ] * df.shape[0]
        if TBcore.get_option('system', 'verbose'):
            sys.stdout.flush()
        return df

    # Get the appropiate path each structure should have.
    data = []
    for _, row in newdf.iterrows():
        data.append(execute(row, structures, flip))
    return pd.concat(data)
    return newdf.apply(execute, axis=1, result_type='broadcast', args=(structures, flip))


def download_pdb( pdbdb: SBIdb.PDBLink,
                  pdb3d: SBIstr.PDBFrame,
                  pdbid: str,
                  pdbchain: str
                  ) -> Tuple[str, SBIstr.PDBFrame]:
    filename = pdbdb.store_local_path('{0}_{1}.pdb.gz'.format(pdbid, pdbchain))
    if not os.path.isfile(filename):
        if TBcore.get_option('system', 'verbose'):
            sys.stdout.write('  Downloading PDB {}\n'.format(pdbid))
        pdb3d = SBIstr.PDB('fetch:{0}'.format(pdbid), format='pdb', clean=True,
                           dehydrate=True, hetatms=False)['Chain:{}'.format(pdbchain)]
        pdb3d.write(filename, format='pdb')
        return filename, pdb3d

    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('  {} already exists\n'.format(filename))
    if pdb3d is None or pdb3d.id != '{0}_{1}'.format(pdbid, pdbchain):
        pdb3d = SBIstr.PDB('fetch:{0}'.format(pdbid), format='pdb', clean=True,
                           dehydrate=True, hetatms=False)['Chain:{}'.format(pdbchain)]
        return filename, pdb3d
    else:
        return filename, pdb3d
