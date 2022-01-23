# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>
.. codeauthor:: Zander Harteveld <zandermilanh@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""

from collections import Iterable
import copy
from random import random
from bisect import bisect
import numpy as np
import scipy.spatial
from transforms3d.euler import euler2mat, mat2euler


class VirtualStructure(object):

    _MAX_AA_DIST = 3.2
    _ATOMTYPE     = ("N", "CA", "C", "O")#, "H")
    _STRING_X    = "HETATM{0:>5d}  X     X {2}{0:>4d} {1[0]:>11.3f}{1[1]:>8.3f}{1[2]:>8.3f}  1.00"

    _STRING_ATOMS   = {
                        "N": "ATOM  {4:>5d}  N   {3:>3} {2}{0:>4d} {1[0]:>11.3f}{1[1]:>8.3f}{1[2]:>8.3f}  1.00",
                        "CA": "ATOM  {4:>5d}  CA  {3:>3} {2}{0:>4d} {1[0]:>11.3f}{1[1]:>8.3f}{1[2]:>8.3f}  1.00",
                        "C": "ATOM  {4:>5d}  C   {3:>3} {2}{0:>4d} {1[0]:>11.3f}{1[1]:>8.3f}{1[2]:>8.3f}  1.00",
                        "O": "ATOM  {4:>5d}  O   {3:>3} {2}{0:>4d} {1[0]:>11.3f}{1[1]:>8.3f}{1[2]:>8.3f}  1.00",
                        "H": "ATOM  {4:>5d}  H   {3:>3} {2}{0:>4d} {1[0]:>11.3f}{1[1]:>8.3f}{1[2]:>8.3f}  1.00",
                        }

    _A123 = {"C": "CYS", "D": "ASP", "S": "SER", "Q": "GLN", "K": "LYS",
             "I": "ILE", "P": "PRO", "T": "THR", "F": "PHE", "N": "ASN",
             "G": "GLY", "H": "HIS", "L": "LEU", "R": "ARG", "W": "TRP",
             "A": "ALA", "V": "VAL", "E": "GLU", "Y": "TYR", "M": "MET"}
    _A321 = {v: k for k, v in _A123.items()}
    _AA_STAT = [("G", 1.00)]
    _TYPE = "C"

    def __init__(self, residues, centre = [0., 0., 0.], chain = "A"):
        self.residues = int(residues)
        self.chain    = chain
        self.centre   = np.array(centre, dtype="float64")
        self.max_dist = float(self._MAX_AA_DIST * self.residues)
        self.edges    = [np.copy(self.centre) + np.array([0, self.max_dist / 2, 0]),
                         np.copy(self.centre) - np.array([0, self.max_dist / 2, 0])]
        self.points   = []
        for x in range(self.residues):
            self.points.append(np.copy(self.edges[0]) - np.array([0, self._MAX_AA_DIST * x, 0]) )
        self.atoms    = []
        self.atomtypes = []
        self.ca_atoms = []
        self.atom = None
        self.Rapplied = np.eye(3)
        self.is_inverted = False
        self.spinned     = 0
        self.sequence    = None
        self.ref         = None
        self.name        = None

    # BOOLEANS
    def in_origin(self):
        return np.allclose(self.centre, [0., 0., 0.])

    def goes_up(self):
        if len(self.atoms) > 0:
            return self.atoms[0][1] < self.atoms[-1][1]
        elif len(self.points) > 0:
            return self.points[0][1] < self.points[-1][1]

    def goes_down(self):
        return not self.goes_up()

    # GETTERS
    def get_type(self):
        return self._TYPE

    def up_is_1(self):
        return 1 if self.goes_up() else 0

    # TILT
    def tilt_x_degrees(self, angle): self.tilt_degrees(x_angle = angle)
    def tilt_y_degrees(self, angle): self.tilt_degrees(y_angle = angle)
    def tilt_z_degrees(self, angle): self.tilt_degrees(z_angle = angle)

    def tilt_degrees(self, x_angle = 0, y_angle = 0, z_angle = 0, store = True):
        if x_angle == 0 and y_angle == 0 and z_angle == 0: return
        self.tilt_radiants(x_angle = np.radians(x_angle),
                           y_angle = np.radians(y_angle),
                           z_angle = np.radians(z_angle), store = store)

    def tilt_x_radiants(self, angle): self.tilt_radiants(x_angle = angle)
    def tilt_y_radiants(self, angle): self.tilt_radiants(y_angle = angle)
    def tilt_z_radiants(self, angle): self.tilt_radiants(z_angle = angle)

    def tilt_radiants(self, x_angle = 0, y_angle = 0, z_angle = 0, store = True):
        Rx = euler2mat(x_angle, 0, 0, "sxyz")
        Ry = euler2mat(0, y_angle, 0, "sxyz")
        Rz = euler2mat(0, 0, z_angle, "sxyz")
        R  = np.dot(Rz, np.dot(Rx, Ry))

        tmpctr = np.array([0., 0., 0.])
        fixpos = not np.allclose(self.centre, tmpctr)
        tmpctr = np.copy(self.centre)
        if fixpos: self.shift(x = -tmpctr[0], y = -tmpctr[1], z = -tmpctr[2])
        self.apply_matrix(R)
        if fixpos: self.shift(x = tmpctr[0], y = tmpctr[1], z = tmpctr[2])

        if store: self.Rapplied = np.dot(self.Rapplied, R)

    def apply_matrix(self, R):
        if len(self.edges):  self.edges  = np.dot(self.edges,  R)
        if len(self.points): self.points = np.dot(self.points, R)
        if len(self.atoms):  self.atoms  = np.dot(self.atoms,  R)

    # ROTATE ON AXIS
    def spin_radians(self, angle): self.spin_degrees(np.degrees(angle))

    def spin_degrees(self, angle):
        if np.allclose(self.Rapplied, np.eye(3)):
            self.tilt_degrees(y_angle = angle, store = False)
        else:
            euler1 = mat2euler(self.Rapplied)
            euler2 = mat2euler(self.Rapplied.transpose())
            self.tilt_radiants(euler2[0], euler2[1], euler2[2])
            self.tilt_degrees(y_angle = angle, store = False)
            self.tilt_radiants(euler1[0], euler1[1], euler1[2])
        self.spinned += angle

    # SHIFT
    def shift_x(self, x): self.shift(x, 0., 0.)
    def shift_y(self, y): self.shift(0., y, 0.)
    def shift_z(self, z): self.shift(0., 0., z)

    def shift(self, x = 0., y = 0., z = 0.):
        t = np.array(x) if isinstance(x, Iterable) else np.array([x, y, z])
        self.centre += t
        self.edges  += t
        self.points += t
        if len(self.atoms): self.atoms += t

    def shift_to_origin(self):
        anti = np.copy(self.centre) if not self.in_origin() else np.array([0., 0., 0.])
        self.shift(-anti)
        return anti

    # FUNC
    def invert_direction(self):
        if np.allclose(self.Rapplied, np.eye(3)):
            self.tilt_degrees(x_angle = 180, y_angle = 180, store = False)
        else:
            euler1 = mat2euler(self.Rapplied)
            euler2 = mat2euler(self.Rapplied.transpose())
            self.tilt_radiants(euler2[0], euler2[1], euler2[2])
            self.tilt_degrees(x_angle = 180, y_angle = 180, store = False)
            self.tilt_radiants(euler1[0], euler1[1], euler1[2])
        self.is_inverted = not self.is_inverted

    def duplicate(self):
        return copy.deepcopy(self)

    def add_3AAseq(self, seq_array):
        self.add_1AAseq([self._A321[x] for x in seq_array])

    def add_1AAseq(self, seq_array):
        if isinstance(seq_array, list):
            self.sequence = "".join(seq_array)
        elif isinstance(seq_array, str):
            self.sequence = seq_array

    def length(self, points = False):
        if not points: return 3
        else:          return len(self.points) + 3

    def create_stat_sequence(self):
        def weighted_choice(choices):
            values, weights = zip(*choices)
            total = 0
            cum_weights = []
            for w in weights:
                total += w
                cum_weights.append(total)
            x = random() * total
            i = bisect(cum_weights, x)
            return values[i]

        if self.sequence is None:
            self.sequence = ""
            for x in range(self.residues):
                self.sequence += weighted_choice(self._AA_STAT)

    def remove_movement_memory(self):
        self.Rapplied = np.eye(3)

    def center_distance_to(self, other):
        return scipy.spatial.distance.euclidean(self.centre, other.centre)

    # PRINT
    def edge_points(self, atom = 1):
        data = []
        data.append(self._STRING_X.format(atom + 0, self.edges[0], self.chain))
        data.append(self._STRING_X.format(atom + 1, self.edges[1], self.chain))
        return "\n".join(data)

    def center_point(self, atom = 1):
        return self._STRING_X.format(atom, self.centre, self.chain)

    def axis_points(self, atom = 1):
        data = []
        data.append(self._STRING_X.format(atom + 0, self.edges[0], self.chain))
        data.append(self._STRING_X.format(atom + 1, self.centre, self.chain))
        data.append(self._STRING_X.format(atom + 2, self.edges[1], self.chain))
        return "\n".join(data)

    def guide_points(self, atom = 1):
        count = 0
        data  = []
        for x in self.points:
            data.append(self._STRING_X.format(atom + count, x, self.chain))
            count += 1
        return "\n".join(data)

    def vector_points(self, atom = 1, points = False):
        data = []
        data.append(self.axis_points(atom))
        if points: data.append(self.guide_points(atom + 3))
        return "\n".join(data)

    def secondary_structure_def(self):
        pass

    def atom_points(self, atom = 1, seq = None):
        count = 0
        data  = []
        d = {"C": "CYS", "D": "ASP", "S": "SER", "Q": "GLN", "K": "LYS",
             "I": "ILE", "P": "PRO", "T": "THR", "F": "PHE", "N": "ASN",
             "G": "GLY", "H": "HIS", "L": "LEU", "R": "ARG", "W": "TRP",
             "A": "ALA", "V": "VAL", "E": "GLU", "Y": "TYR", "M": "MET"}
        if seq is None and self.sequence is not None: seq = self.sequence #* (int(len(self.atoms)/len(seq)))
        if seq is None: seq = "G" #* (int(len(self.atoms)/len(seq)))
        else:           seq = seq.upper() #* (int(len(self.atoms)/len(seq)))

        for x, (points, atomtype) in enumerate(zip(self.atoms, self.atomtypes)):
            data.append(self._STRING_ATOMS[atomtype].format(atom + count, points, self.chain, d[seq[count]], atom + x))
            if (1 + x)%len(self._ATOMTYPE)==0:
                count += 1
        return "\n".join(data)

    def __eq__(self, other):
        # TODO this should be better...
        return self.centre == other.centre

    def __ne__( self, other ):
        return not self.__eq__( other )

    def __hash__( self ):
        return str(self.centre)

if __name__ == '__main__':
    b = VirtualStructure(15, [0, 0, 0], chain="Y")
    print (b.axis_points(1))
    x = VirtualStructure(15, [0, 0, 0], chain="X")
    x.tilt_degrees(x_angle = 90)
    print (x.axis_points(4))
    z = VirtualStructure(15, [0, 0, 0], chain="Z")
    z.tilt_degrees(z_angle = 90)
    print (z.axis_points(7))
    # AXIS
    p1 = VirtualStructure(15, [0, 0, 0], chain="A")
    p1.tilt_degrees(z_angle=45)
    print (p1.axis_points(10))
    p2 = VirtualStructure(15, [0, 0, 0], chain="B")
    p2.tilt_degrees(x_angle=45)
    print (p2.axis_points(13))
    p3 = VirtualStructure(15, [0, 0, 0], chain="C")
    p3.tilt_degrees(x_angle=45, z_angle=45)
    print (p3.axis_points(16))
    p3 = VirtualStructure(15, [0, 0, 0], chain="D")
    p3.tilt_degrees(x_angle=45, y_angle=90, z_angle=45)
    print (p3.axis_points(19))
    # TURN ON CENTER MASS

    b.chain = "E"
    b.shift(3., 12., 6.)
    print (b.axis_points(22))
    x = VirtualStructure(15, [3., 12., 6.], chain="F")
    x.tilt_degrees(x_angle = 90)
    print (x.axis_points(25))
    z = VirtualStructure(15, [3., 12., 6.], chain="G")
    z.tilt_degrees(z_angle = 90)
    z.shift(x=4.9, y=1., z=0.)
    print (z.axis_points(28))

    y = VirtualStructure(15, [3., 12., 6.], chain="C")
    x = VirtualStructure(15, [3., 12., 6.], chain="C")
    x.tilt_degrees(z_angle=45)
    z = VirtualStructure(15, [3., 12., 6.], chain="D")
    z.tilt_degrees(x_angle=45)
    p = VirtualStructure(15, [3., 12., 6.], chain="E")
    p.tilt_degrees(x_angle=45, z_angle=45)
    p.tilt_x_degrees(45)
    p.tilt_y_degrees(45)
    print y.axis_points(19)
    print x.axis_points(22)
    print z.axis_points(25)
    print p.axis_points(28)
    t = VirtualStructure(15, [3., 12., 6.], chain="C")
    t.tilt_degrees(x_angle = 45, z_angle = 90)
    t.invert_direction()
    print t.axis_points(31)
