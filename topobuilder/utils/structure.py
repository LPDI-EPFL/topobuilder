# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>
.. codeauthor:: Zander Harteveld <zandermilanh@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import List, Dict, Tuple, Union, Optional
from pathlib import Path
import sys
import math
from string import ascii_uppercase

# External Libraries
try:
    from SBI.structure import PDB, Frame3D
    import SBI.core as SBIcr
except ImportError:
    class Frame3D():
        pass

from logbook import Logger
import pandas as pd
import numpy as np
import sympy as sy

# This Library
import topobuilder.core as TBcore


__all__ = ['build_pdb_object', 'get_loop_length', 'pdb_geometry_from_rules', 'reverse_motif', 'StructuralError']


np.set_printoptions(precision=3)


def build_pdb_object( log: Logger,
                      sses: List[Dict],
                      loops: Union[List[int], int],
                      concat: Optional[bool] = True,
                      outfile: Optional[Union[str, Path]] = None
                      ) -> Tuple[Frame3D, List[int]]:
    """Make the parametrically build atoms in a :class:`.Case` into a PDB file.

    :param log: Job logger.
    :param sses: List of the secondary structures to build. Each SSE dictionary must contain the
        ``metadata.atoms`` keys, already in the final expected position.
    :param loops: Number of residues between SSE. It can be one less than the number of structures,
        which assumes no N- or C-terminal, or one more, which assumes N- and C-terminal residues.
    :param concat: When :data:`True`, return the full stucture as a single object, otherwise
        return a list of the individual parts.
    :param outfile: If provided, write the structure to file.
    """
    if isinstance(loops, int):
        loops = [loops, ] * (len(sses) - 1)

    if len(loops) != len(sses) - 1:
        raise ValueError('Number of loops should equal number of SSE minus one.')

    pieces = []
    columns = ['auth_comp_id', 'auth_atom_id', 'auth_seq_id', 'Cartn_x', 'Cartn_y', 'Cartn_z']
    start = 1 if len(loops) < len(sses) else loops.pop(0)
    log.debug(f'starting numbering with: {start}')
    for i, sse in enumerate(sses):
        start = start if i == 0 else int(sses[i - 1]['length']) + loops[i - 1] + start
        pdb_numbering = pd.DataFrame(sse['metadata']['atoms'], columns=columns)['auth_seq_id'].values
        try:
            structure = PDB(pd.DataFrame(sse['metadata']['atoms'], columns=columns)).renumber(start)
        except:
            structure = PDB(pd.DataFrame(sse['metadata']['atoms'], columns=columns))
            structure['auth_seq_id'] += (start - structure['auth_seq_id'].values[0])

        structure = structure.assign(sse_id=[sse["id"]]*len(structure), pdb_num=pdb_numbering)
        pieces.append(structure)

    structure = pd.concat(pieces, sort=False).reset_index()
    structure['id'] = list(range(1, structure.shape[0] + 1))

    if outfile is not None:
        structure.write(output_file=str(outfile), format='pdb', clean=True,
                        force=TBcore.get_option('system', 'overwrite'))

    if not concat:
        return pieces

    return structure, [int(p.iloc[-1]['auth_seq_id']) for p in pieces]


def get_loop_length( log: Logger, sse1: Frame3D, sse2: Frame3D, loop_step: int, loop_range: int) -> Tuple[int, int]:
    """Calculate the expected number of residues to join two SSE.

    :param log: Job Logger.
    :param sse1: N-SSE.
    :param sse2: C-SSE.
    :param loop_step: Assumption on how much distance a residue can cover.
    :param loop_range: Plus-minus range of residue length.
    """
    from SBI.structure import ChainFrame
    from SBI.structure.geometry.basics import distance

    res1 = ChainFrame(PDB(sse1)).last_compound
    res2 = ChainFrame(PDB(sse2)).first_compound
    distance = distance(res1[res1['label_atom_id'] == 'C'].coordinates,
                        res2[res2['label_atom_id'] == 'N'].coordinates)
    log.debug(f'Distance between SSE is {distance} Angstrongs.')
    distance = math.ceil(distance / loop_step)
    log.debug(f'Assuming the need of {distance} residues with a {loop_range} residue range.')
    distance = [x for x in range(distance - loop_range - 1, distance + loop_range + 1) if x > 0]
    return max(distance), min(distance)


def pdb_geometry_from_rules( pdb_file: Union[Path, str, Frame3D],
                             rules: List[Tuple],
                             log: Optional[Logger] = None
                            ) -> pd.DataFrame:
    """Calculates the geometry statistic from a PDB.

    :param log: Job Logger.
    :param pdb_file: The pdb file to calculate the geometry from.
    :param rules: The rules to be applied.
    """
    if isinstance(pdb_file, (Path, str)):
        pdb_file = Path(pdb_file)
        if not pdb_file.is_file():
            raise IOError('PDB structure {} not found.'.format(pdb_file))
        pdb3d = PDB(str(pdb_file), format='pdb', clean=True, dehydrate=True, hetatms=False)['AtomTask:PROTEINBACKBONE']
    elif isinstance(pdb_file, Frame3D):
        pdb3d = pdb_file['AtomTask:PROTEINBACKBONE']
    else:
        raise ValueError('Unexpected type for pdb_file.')

    if log:
        log.info(f'PDB:Analyzing geometry of {pdb3d.id}\n')
        log.debug(f'PDB:Available secondary structures {",".join([x[0] for x in rules])}\n')
        log.debug(f'PDB:With ranges {",".join(["{}-{}".format(*x[1]) for x in rules])}\n')
        log.debug(f'PDB:With flip policy {",".join([str(x[2]) for x in rules])}\n')
    else:
        sys.stdout.write(f'PDB:Analyzing geometry of {pdb3d.id}\n')
        sys.stdout.write(f'PDB:Available secondary structures {",".join([x[0] for x in rules])}\n')
        sys.stdout.write(f'PDB:With ranges {",".join(["{}-{}".format(*x[1]) for x in rules])}\n')
        sys.stdout.write(f'PDB:With flip policy {",".join([str(x[2]) for x in rules])}\n')

    pieces = make_pieces(pdb3d, rules)
    pieces = make_vectors(pieces, rules)
    pieces = make_planes(pieces)
    df = make_angles_and_distances(pieces)
    df = df.assign(pdb_path=[str(pdb_file) if not isinstance(pdb_file, Frame3D) else pdb3d.id, ] * df.shape[0])
    return df


def make_pieces( pdb3d: Frame3D, rules: List[Tuple] ) -> Dict:
    """Chunks a PDB into its SSE pieces.

    :param pdb3d: PDB stored as a frame.
    :param rules: The rules to be applied.
    """
    pieces = {}
    for piece in rules:
        sse_id, ranges, _ = piece
        if TBcore.get_option('system', 'debug'):
            sys.stdout.write('PDB:Individualize SSE:{} of {}\n'.format(sse_id, pdb3d.id))
        with SBIcr.on_option_value('structure', 'source', 'label'):
            segment = pdb3d['Residue:{0}-{1}'.format(ranges[0], ranges[1])]
            pieces.setdefault(sse_id, {}).setdefault('atoms', segment)
        with SBIcr.on_option_value('structure', 'source', 'auth'):
            if TBcore.get_option('system', 'verbose'):
                try:
                    first, last = segment.first_compound.number, segment.last_compound.number
                    sys.stdout.write('PDB:{2} - Range: {0}-{1}\n'.format(first, last, pdb3d.id))
                except:
                    continue
    return pieces


def make_vectors( pieces: Dict, rules: List[Tuple] ) -> Dict:
    """Calculates the vectors for a SSE pieces.

    :param pieces: The SSE pieces to calculate the vectors from.
    :param rules: The rules to be applied.
    """
    for piece in rules:
        sse_id, _, flip = piece
        if TBcore.get_option('system', 'debug'):
            sys.stdout.write('PDB:Getting vector for SSE:{}{}\n'.format(sse_id, '' if not flip else ' - flipped!'))
        pieces[sse_id].setdefault('vector', list(pieces[sse_id]['atoms'].eigenvectors(40)[-1]))
        if flip:
            pieces[sse_id]['vector'] = [list(x) for x in np.flip(np.asarray(pieces[sse_id]['vector']), axis=0)]
        print('PYMOL:make_arrow("vect{0}", {1}, {2})'.format(sse_id, np.asarray(pieces[sse_id]['vector'][0]).tolist(),
                                                             np.asarray(pieces[sse_id]['vector'][-1]).tolist()))
    return pieces


def make_planes( pieces: Dict ) -> Dict:
    """Calculates the planes for SSE pieces.

    :param pieces: The SSE pieces to calculate the vectors from.
    """
    blayers = sorted(set([x[0] for x in pieces if x.endswith('E') and len(x) == 3]))
    hlayers = [x[0] for x in pieces if x.endswith('H') and len(x) == 3]
    hlayers = sorted(set([x for x in hlayers if hlayers.count(x) > 1]))

    for layer in blayers:
        if TBcore.get_option('system', 'debug'):
            sys.stdout.write('PDB:Generating plane for beta layer {}\n'.format(layer))
        pieces.setdefault(layer, {'layer': [], 'floor': [], 'side': []})
        structure = pd.concat([pieces[x]['atoms'] for x in pieces if len(x) == 3 and x.startswith(layer)])
        eign = structure.eigenvectors(30)
        eign = eigenlayers_fix(eign, [pieces[x]['vector'] for x in pieces if len(x) == 3 and x.startswith(layer)])
        pieces[layer]['layer'] = [list(eign[1][0]), list(eign[1][-1]), list(eign[2][0])]
        pieces[layer]['floor'] = [list(eign[1][0]), list(eign[1][-1]), list(eign[0][0])]
        pieces[layer]['side'] = [list(eign[2][0]), list(eign[2][-1]), list(eign[0][0])]
        for i, v in enumerate([('floor', 'red'), ('side', 'green'), ('layer', 'blue')]):
            print('PYMOL:make_arrow("b{0}", {1}, {2}, color="{3}")'.format(v[0], eign[i][0].tolist(),
                                                                           eign[i][-1].tolist(), v[1]))

    for layer in hlayers:
        if TBcore.get_option('system', 'debug'):
            sys.stdout.write('PDB:Generating plane for helix layer {}\n'.format(layer))
        pieces.setdefault(layer, {'layer': [], 'floor': [], 'side': []})
        structure = pd.concat([pieces[x]['atoms'] for x in pieces if len(x) == 3 and x.startswith(layer)])
        eign = structure.eigenvectors(30)
        eign = eigenlayers_fix(eign, [pieces[x]['vector'] for x in pieces if len(x) == 3 and x.startswith(layer)])
        pieces[layer]['layer'] = [list(eign[1][0]), list(eign[1][-1]), list(eign[2][0])]
        pieces[layer]['floor'] = [list(eign[1][0]), list(eign[1][-1]), list(eign[0][0])]
        pieces[layer]['side'] = [list(eign[2][0]), list(eign[2][-1]), list(eign[0][0])]
        for i, v in enumerate([('floor', 'red'), ('side', 'green'), ('layer', 'blue')]):
            print('PYMOL:make_arrow("h{0}", {1}, {2}, color="{3}")'.format(v[0], eign[i][0].tolist(),
                                                                           eign[i][-1].tolist(), v[1]))

    return pieces


def eigenlayers_fix( eign: np.ndarray, vectors: np.ndarray, scape: bool = False ) -> np.ndarray:
    """Orients the plane calculate from the pieces.
    """
    eign = eign.copy()
    eign2 = sy.Line(eign[2][0], eign[2][-1])
    angles_layer = []
    for sse in vectors:
        angles_layer.append(math.degrees(sy.Line(sse[0], sse[-1]).angle_between(eign2)))

    # Fix layer direction
    if len(vectors) > 3:
        do_flip = np.isclose(angles_layer, [180, ] * len(angles_layer), atol=40)
        do_flip = sum(do_flip) >= len(do_flip) - 1
    else:
        do_flip = np.allclose(angles_layer, [180, ] * len(angles_layer), atol=40)
    # Fix vectors -> switch side and layer
    if len(vectors) > 3:
        do_switch = [not x for x in np.isclose(angles_layer, [0, ] * len(angles_layer), atol=40)]
        do_switch = sum(do_switch) >= len(do_switch) - 1
    else:
        do_switch = not np.allclose(angles_layer, [0, ] * len(angles_layer), atol=35)

    # Apply
    if do_flip:
        eign[2] = [list(x) for x in np.flip(np.asarray(eign[2]), axis=0)]
    elif do_switch and not scape:
        eign[[1, 2]] = eign[[2, 1]]
        eign = eigenlayers_fix(eign, vectors, True)

    return eign


def make_angles_and_distances( pieces: Dict ) -> pd.DataFrame:
    """Calculates the angles and distances from the vectors and planes.

    :param pieces: The SSE pieces to calculate the vectors from.
    """
    data = {'sse': [], 'layer': [],
            'angles_layer': [], 'angles_floor': [], 'angles_side': [],
            'points_layer': [], 'points_floor': [], 'points_side': [],
            'tilted_layer': [], 'tilted_floor': [], 'tilted_side': []}

    for layer in sorted(set([x[0] for x in pieces if len(x) == 1])):
        for sse in [x for x in pieces if len(x) == 3]:
            if abs(ascii_uppercase.find(layer) - ascii_uppercase.find(sse[0])) <= 1:
                data['sse'].append(sse)
                data['layer'].append(layer)
                for iplane, plane in enumerate(pieces[layer]):
                    if TBcore.get_option('system', 'debug'):
                        sys.stdout.write('PDB:{} geometry plane {} vs. sse {}\n'.format(plane, layer, sse))
                    syPlane = sy.Plane(sy.Point3D(pieces[layer][plane][0]),
                                       sy.Point3D(pieces[layer][plane][1]),
                                       sy.Point3D(pieces[layer][plane][2]))
                    syLine = sy.Line(pieces[sse]['vector'][0], pieces[sse]['vector'][-1])
                    syPoint = sy.Point3D(*pieces[sse]['vector'][1])
                    data[f'angles_{plane}'].append(math.degrees(syPlane.angle_between(syLine)))
                    data[f'points_{plane}'].append(float(syPlane.distance(syPoint)))
                    data[f'tilted_{plane}'].append(float(syPlane.distance(default_plane(iplane))))
    return pd.DataFrame(data)


def default_plane( pick: int ) -> sy.Plane:
    """
    """
    x = sy.Point3D([30, 0, 0])
    y = sy.Point3D([0, 30, 0])
    z = sy.Point3D([0, 0, 30])
    c = sy.Point3D([0, 0, 0])

    if pick == 0:
        return sy.Plane(y, c, x)
    elif pick == 1:
        return sy.Plane(x, c, z)
    elif pick == 2:
        return sy.Plane(z, c, y)
    else:
        raise ValueError('Selection must be between 0 and 2')


def reverse_motif( log: Logger,
                   source: Union[str, Path],
                   selection: List[str],
                   attach: List[str],
                   hotspot: str,
                   identifier: str,
                   binder: Optional[str] = None, ) -> Dict:
    """Process a provided motif so that it can be attached to a :term:`FORM`.

    :param log: Logger from the calling :class:`.Node` to keep requested verbosity.
    :param source: File containing the structural data.
    :param selection: Selection defining the motif of interest.
    :param hotspot: Single position defining the exposed side.
    :param identifier: Single position defining the exposed side.
    :param binder: Selection defining the binder.

    """
    # Load Structure
    pdbSTR = PDB(source, header=False, dehydrate=True)
    # Get the full motif to define its planes
    motif, hotspots  = pick_motif(log, pdbSTR, selection, attach, hotspot)
    # Pick Binder
    binder = pick_binder(log, pdbSTR, binder)

    # # Find motif's orientation
    # eigens = dict(map(reversed, zip(motif['AtomType:CA'].eigenvectors(10), ('perpendicular', 'side', 'major'))))
    # edist = [np.linalg.norm(hotspot.coordinates - eigens['perpendicular'][0]), np.linalg.norm(hotspot.coordinates - eigens['perpendicular'][-1]),
    #          np.linalg.norm(hotspot.coordinates - eigens['side'][0]), np.linalg.norm(hotspot.coordinates - eigens['side'][-1])]
    # mdist = edist.index(min(edist))
    # if mdist > 1:  # swap side and perpendicular
    #     tmp = eigens['side']
    #     eigens['side'] = eigens['perpendicular']
    #     eigens['perpendicular'] = tmp
    #     if mdist == 2:  # Change perpendicular orientation
    #         eigens['perpendicular'] = np.flip(eigens['perpendicular'], axis=0)
    # if mdist == 0:  # Change perpendicular orientation
    #     eigens['perpendicular'] = np.flip(eigens['perpendicular'], axis=0)
    #
    # # Try to identify the orientation of the motif.
    # for k in eigens:
    #     pymol_arrow(f'arrow_{k}', eigens[k][0], eigens[k][-1],
    #                 'white' if k == 'major' else 'red' if k == 'side' else 'green')
    return motif, binder, hotspots, attach, selection, identifier


def pick_motif( log: Logger, pdbSTR: Frame3D, selection: List[str], attach: List[str], hotspot: str ) -> Dict:
    """Saves the motif from a given PDB.

    :param log: Logger from the calling :class:`.Node` to keep requested verbosity.
    :param pdbSTR: The PDB to take the motif from.
    :param selection: Residues to be considered as motif.
    :param attach: The :term:`SKETCH` part to be employed to add the motif onto.
    :param hotspot: The hotspot residues on the motif to keep fixed during design.
    """
    # Get the motif segments
    log.info(f'Loading motifs at {",".join(selection)}.')
    segments = []
    identification = []

    for s, a in zip(selection, attach):
        s = s.split(':')
        st = pdbSTR[f'Chain:{s[0]}'][f'Residue:{s[1]}']
        segments.append(st)
        identification.extend([a] * len(st))
    motif = pd.concat(segments)

    motif = motif.assign(sse_id=identification,
                         internal_num=motif['auth_seq_id'] - (motif['auth_seq_id'].values[0] - 1))

    # Get the non-segments
    # log.info('Loading regions between segments.')
    # midsegs = []
    # for c, i in enumerate(range(1, len(shape), 2)):
    #     log.debug(f'Current separator is {shape[i]}')
    #     log.debug(f'Placed between selections {selection[c]} and {selection[c + 1]}')
    #     if shape[i] == 'x':
    #         midsegs.append(None)
    #     else:
    #         chain = [selection[c].split(':')[0], selection[c + 1].split(':')[0]]
    #         if chain[0] != chain[1]:
    #             raise StructuralError('Continuos motif cannot be picked from different chains.')
    #         resi = [selection[c].split(":")[-1].split("-")[-1], selection[c + 1].split(":")[-1].split("-")[0]]
    #         midsegs.append(pdbSTR[f'Chain:{chain[0]}'][f'Residue:{resi[0]}-{resi[1]}'])
    #         midsegs[-1] = midsegs[-1].pop_out(0).pop_out(-1)

    # Get the orientation-guiding hotspot
    log.info(f'Picking the hotspot residues {hotspot}.')
    htspot = pick_hotspots(log, motif, hotspot)

    motif = list(motif[['auth_comp_id', 'auth_atom_id', 'auth_seq_id', 'auth_asym_id', 'sse_id', 'internal_num',
                        'Cartn_x', 'Cartn_y', 'Cartn_z']].values)
    # htspot = motif[f'Chain:{hotspot.split(":")[0]}'][f'Residue:{hotspot.split(":")[1]}']['AtomType:CA']
    # if htspot.is_empty:
    #     raise StructuralError('Orientation-guiding hotspot not found inside the motif.')
    # if len(htspot) > 1:
    #     raise StructuralError('Orientation-guiding hotspot should be a single residue.')

    # Get the orientation
    #eigens = single_segment_eigens(log, motif, htspot) if len(segments) == 1 else multi_segment_eigens(log, segments, hotspot)

    # Try to identify the orientation of the motif.
    # for i, e in enumerate(eigens):
    #     for k in e:
    #         pymol_arrow(f'arrow_{motif.id}_{i + 1}_{k}', e[k][0], e[k][-1],
    #                     'white' if k == 'major' else 'red' if k == 'side' else 'green')

    return motif, htspot


def single_segment_eigens( log: Logger, motif: Frame3D, hotspot: Frame3D ) -> List[Dict]:
    """Generate eigenvectors for a single segment motif.

    :param motif: motif to geometrically analyze.
    :param hotspot: CA atom to guide front-view.
    """
    log.debug('Single segment geometric analysis.')
    # Get the orientation
    eigens = dict(map(reversed, zip(motif['AtomType:CA'].eigenvectors(10), ('perpendicular', 'side', 'major'))))
    edist = [np.linalg.norm(hotspot.coordinates - eigens['perpendicular'][0]), np.linalg.norm(hotspot.coordinates - eigens['perpendicular'][-1]),
             np.linalg.norm(hotspot.coordinates - eigens['side'][0]), np.linalg.norm(hotspot.coordinates - eigens['side'][-1])]
    mdist = edist.index(min(edist))
    if mdist > 1:  # swap side and perpendicular
        log.debug('Swap side and perpendicular axes.')
        tmp = eigens['side']
        eigens['side'] = eigens['perpendicular']
        eigens['perpendicular'] = tmp
        if mdist == 2:  # Change perpendicular orientation
            log.debug('Flip perpendicular axis.')
            eigens['perpendicular'] = np.flip(eigens['perpendicular'], axis=0)
    if mdist == 0:  # Change perpendicular orientation
        log.debug('Flip perpendicular axis.')
        eigens['perpendicular'] = np.flip(eigens['perpendicular'], axis=0)
    return [eigens, ]


def multi_segment_eigens( log: Logger, motifs: List[Frame3D], hotspot: str ) -> List[Dict]:
    """Generate eigenvectors for a multisegment motif.

    :param motifs: list of segments of the motif to geometrically analyze.
    :param hotspot: CA atom to guide front-view.
    """
    log.debug('Multiple segment geometric analysis.')
    eigens = []
    hotspot_segment = 0
    for i, motif in enumerate(motifs):
        eigens.append(dict(map(reversed, zip(motif['AtomType:CA'].eigenvectors(10), ('perpendicular', 'side', 'major')))))
        if not motif[f'Chain:{hotspot.split(":")[0]}'][f'Residue:{hotspot.split(":")[1]}'].is_empty:
            log.debug(f'Guiding residue located in segment {i + 1}.')
            hotspot_segment = i
    return eigens


def pick_binder( log: Logger, pdbSTR: Frame3D, binder: Optional[str] = None ) -> Optional[Frame3D]:
    """Pick binder selection.
    """
    # Pick binder if any
    if binder is None:
        return None
    bndr = []
    for b in binder.split(','):
        if len(b) == 1:
            bndr.append(pdbSTR[f'Chain:{b}'])
        else:
            b = b.split(':')
            bndr.append(pdbSTR[f'Chain:{b[0]}'][f'Residue:{b[1]}'])
    bndr = pd.concat(bndr)
    bndr = list(bndr[['auth_comp_id', 'auth_atom_id', 'auth_seq_id', 'auth_asym_id',
                      'Cartn_x', 'Cartn_y', 'Cartn_z']].values)
    return bndr


def pick_hotspots( log: Logger, pdbSTR: Frame3D, hotspots: Optional[str] = None ) -> Optional[Frame3D]:
    """Pick binder selection.
    """
    # Pick binder if any
    if hotspots is None:
        return None
    htspts = []
    for h in hotspots.split(','):
        h = h.split(':')
        htspts.append(pdbSTR[f'Chain:{h[0]}'][f'Residue:{h[1]}'])
    return pd.concat(htspts)


def pymol_arrow( name: str, ini: np.ndarray, end: np.ndarray, color='white' ) -> None:
    """
    """
    print(f'make_arrow {name}, {list(ini)}, {list(end)}, color={color}')


def pymol_point( name: str, coords: np.ndarray ) -> None:
    """
    """
    print(f'pseudoatom {name}, pos={list(coords)}')


class StructuralError(Exception):
    """
    """
