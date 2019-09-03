# -*- coding: utf-8 -*-
# @Author: bonet
# @Date:   2016-04-28 13:49:48
# @Last modified by:   bonet
# @Last modified time: 03-Sep-2019
import json
import numpy as np


def load_json(filename):
    with open(filename) as fd:
        data = json.loads(''.join([l.strip() for l in fd]))
    if "config" not in data: raise AttributeError("config info is mandatory")
    if "name" not in data["config"]: raise AttributeError("config 'name' is mandatory")
    if "status" not in data["config"]: data["config"]["status"] = 0
    if "vall" not in data["config"]: raise AttributeError("provide the path to the vall database to create the fragments")
    if "rbin" not in data["config"]: raise AttributeError("provide the path to the rosetta binaries as 'rbin'")
    return data


def print_json(data, filename = None):
    d = json.dumps(data, indent=2, separators=(',', ': '))
    if filename is None: print d
    else:
        with open(filename, "w") as fd: fd.write(d)


def read_pdbs(data):
    def good_chain(line, chain):
        if not line.startswith("ATOM"): return False
        if line[13:15].strip() not in ["CA", "O", "N", "C"]:        return False
        #if line[13:15].strip() not in ["CA"]:        return False
        if line[21] != chain:           return False
        return True

    def good_line(line, chain, res):
        if not good_chain(line, chain):       return False
        if int(line[22:26]) not in res: return False
        return True

    if "motifs" not in data: return
    if data["config"]["status"] >= 1: return

    all_atoms = []
    for cnt, mtf in enumerate(data["motifs"]):
        if "pdbfile" not in mtf:  raise AttributeError("Motifs need a PDB file assigned")
        if "chain" not in mtf:    mtf["chain"] = "A"
        if "segments" not in mtf: raise AttributeError("Segments of the motif need to be specified")
        if "id" not in mtf:       mtf["id"] = "motif" + str(cnt + 1)
        if "lookZ" not in mtf:    mtf["lookZ"] = 1

        segments = set()
        for sgm in mtf["segments"]:
            for i in range(int(sgm["ini"]), int(sgm["end"]) + 1):
                segments.add(i)

        with open(mtf["pdbfile"]) as fd:
            for line in fd:
                if good_chain(line, mtf["chain"]):
                    all_atoms.append(np.array([float(line[27:38].strip()),
                                               float(line[38:46].strip()),
                                               float(line[46:54].strip())],
                                     dtype = np.float64))
                if good_line(line, mtf["chain"], segments):
                    for i, sgm in enumerate(mtf["segments"]):
                        if int(line[22:26]) >= int(sgm["ini"]):
                            if int(line[22:26]) <= int(sgm["end"]):
                                x = i
                    mtf["segments"][x].setdefault("sequence", []).append(line[17:20])
                    mtf["segments"][x].setdefault("coordinates", [])
                    mtf["segments"][x]["coordinates"].append((line[27:38].strip(),
                                                              line[38:46].strip(),
                                                              line[46:54].strip()))
            mtf["core"] = np.mean(all_atoms, axis = 0, dtype = np.float64).tolist()

    data["config"]["status"] = 1


def read_all_atom_pdb(data):
    def good_chain(line, chain):
        if not line.startswith("ATOM"): return False
        if line[21] not in chain:       return False
        return True

    def good_line(line, chain, res = None):
        if not good_chain(line, chain): return False
        if res is not None and int(line[22:26]) not in res: return False
        return True

    class Section(object):
        """docstring for Section"""
        def __init__(self, ref):
            self.ref = ref
            self.guide = []
            self.atoms = []
            self.amino = []
            self.binat = {}
            self.binaa = {}

        def process(self):
            self.guide = np.array(self.guide, dtype = np.float64)
            self.atoms = np.array(self.atoms, dtype = np.float64)
            for k in self.binat:
                self.binat[k] = np.array(self.binat[k], dtype = np.float64)

        def align_to(self, atoms):
            self.atoms = align_by_Kabsch(self.atoms, atoms, self.guide)
            if len(self.binaa):
                self.binat = align_by_Kabsch(self.binat, atoms, self.guide)
            self.to_pdb()

        def to_pdb(self, atom = 1, residue = 1):
            _STRING_CA   = "ATOM  {0:>5d}{3[0]} {3[1]:>3} {2}{4:>4d} {1[0]:>11.3f}{1[1]:>8.3f}{1[2]:>8.3f}  1.00"
            out = []
            for cx, x in enumerate(self.atoms):
                n = self.amino[cx][2] - self.amino[0][2] + residue
                out.append(_STRING_CA.format(atom + cx, x, "A", self.amino[cx], n))
            return "\n".join(out)

        def bind_to_pdb(self, atom = 1, bind_chains = "B"):
            _STRING_CA   = "ATOM  {0:>5d}{3[0]} {3[1]:>3} {2}{4:>4d} {1[0]:>11.3f}{1[1]:>8.3f}{1[2]:>8.3f}  1.00"
            out = []
            for k in self.binat:
                for cx, x in enumerate(self.binat[k]):
                    n = self.binaa[cx][2] - self.binaa[k][0][2] + 1
                    out.append(_STRING_CA.format(atom + cx, x, bind_chains, self.binaa[k][cx], n))
                out.append("TER")
                bind_chains = chr(ord(bind_chains) + 1)
            return "\n".join(out)

    content = []
    for cnt, mtf in enumerate(data["motifs"]):
        for sgm in mtf["segments"]:
            segments = set()
            for i in range(int(sgm["ini"]), int(sgm["end"]) + 1):
                segments.add(i)
            s = Section(mtf["id"] + "." + sgm["id"])
            with open(mtf["pdbfile"]) as fd:
                for line in fd:
                    if good_line(line, [mtf["chain"]], segments):
                        s.amino.append((line[11:16], line[17:20], int(line[22:26])))
                        x, y, z = float(line[27:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())
                        s.atoms.append([x, y, z])
                        if line[13:15] == "CA": s.guide.append([x, y, z])
                    elif "binders" in mtf and good_line(line, mtf["binders"].split()):
                        s.binaa.setdefault(line[21], [])
                        s.binat.setdefault(line[21], [])
                        s.binaa[line[21]].append((line[11:16], line[17:20], int(line[22:26])))
                        x, y, z = float(line[27:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())
                        s.binat[line[21]].append(np.array([x, y, z], dtype = np.float64))
            s.process()
            content.append(s)
    return content


def align_by_Kabsch(moving_selection, target_guide, moving_guide = None):

    if moving_guide is None: moving_guide = np.copy(moving_selection)
    assert len(target_guide) == len(moving_guide)

    moving_center = np.mean(moving_guide, axis = 0, dtype = np.float64)
    target_center = np.mean(target_guide, axis = 0, dtype = np.float64)

    M = kabsch(moving_guide - moving_center, np.copy(target_guide) - target_center)

    moving_selection -= moving_center
    moving_selection = np.dot(moving_selection, M)
    moving_selection += target_center

    return moving_selection


def kabsch(P, Q):
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
