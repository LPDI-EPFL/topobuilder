# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>
.. codeauthor:: Zander Harteveld <zandermilanh@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""

import numpy as np
from VirtualStructure import VirtualStructure as VS


class VirtualHelix(object):
    def __new__(self, *args, **kwargs):
        if "type" in kwargs:
            hid = kwargs["type"].lower()
            del(kwargs["type"])
            if hid in ["alpha", "h"]: return VirtualHelixAlpha(*args, **kwargs)
            if hid in ["310", "g"]:   return VirtualHelix310(*args, **kwargs)
            if hid in ["pi", "i"]:    return VirtualHelixPi(*args, **kwargs)
        else: return VirtualHelixAlpha(*args, **kwargs)


class VirtualHelixAlpha(VS):

    _MAX_AA_DIST  = 1.5
    _ATOMTYPES    = ("N", "CA", "C", "O")#, "H")

    _ATOM_CA_DIST = {"N": 0.841, "CA": 0, "C": -1.029, "O": -2.248, "H": 1.839}
    _RADIUS       = {"N": 1.5, "CA": 2.3, "C": 1.8, "O": 2.1, "H": 1.5}
    _ANGLES       = {"N": -28.3, "CA": 100, "C": 28.9, "O": 24.5, "H": -22.5}

    # CHOP780201 alpha-helix propensity AAindex (Chou-Fasman, 1978b)
    # TO 0: G -> 0.57; P -> 0.57; C -> 0.70
    _AA_STAT = [("A", 1.42), ("L", 1.21), ("R", 0.98), ("K", 1.16), ("N", 0.67),
                ("M", 1.45), ("D", 1.01), ("F", 1.13), ("C", 0.00), ("P", 0.00),
                ("Q", 1.11), ("S", 0.77), ("E", 1.51), ("T", 0.83), ("G", 0.00),
                ("W", 1.08), ("H", 1.00), ("Y", 0.69), ("I", 1.08), ("V", 1.06)]

    _TYPE = "H"

    def __init__(self, residues, centre = [0., 0., 0.], chain = "A"):
        super(VirtualHelixAlpha, self).__init__(residues, centre, chain)
        self.edge_angles = [0., 0.]
        self.atoms = []
        self.atomtypes = []
        self.residuenumbers = []
        self.ca_atoms = []
        count = 0
        for x in range(len(self.points)):
            count += 1
            for atomtype in self._ATOMTYPES:
                if atomtype is "CA":
                    angle = self._ANGLES["CA"] * x
                else:
                    angle = self._ANGLES["CA"] * x + self._ANGLES[atomtype]
                point = np.copy(self.points[x]) + np.array([self._RADIUS[atomtype], self._ATOM_CA_DIST[atomtype], 0.])
                self._tilt_y_point_from_centre(self.points[x], point, np.radians(angle))
                self.atomtypes.append(atomtype)
                self.atoms.append(point)
                self.residuenumbers.append(count)
                if atomtype is "CA":
                    self.ca_atoms.append(point)

    def grow_nterm(self, residues):   raise NotImplementedError
    def shrink_nterm(self, residues): raise NotImplementedError
    def grow_cterm(self, residues):   raise NotImplementedError
    def shrink_cterm(self, residues): raise NotImplementedError

    def secondary_structure_def(self):
        pass

    def _tilt_y_point_from_centre(self, centre, point, angle):
        tmp_point = point[0] - centre[0] , point[2] - centre[2]
        tmp_point = ( tmp_point[0] * np.cos(angle) - tmp_point[1] * np.sin(angle),
                      tmp_point[1] * np.cos(angle) + tmp_point[0] * np.sin(angle))
        tmp_point = tmp_point[0] + centre[0] , tmp_point[1] + centre[2]
        point[0]  = tmp_point[0]
        point[2]  = tmp_point[1]


class VirtualHelix310(VirtualHelixAlpha):

    _MAX_AA_DIST = 2.0
    _RADIUS      = 1.9
    _TYPE        = "G"

    def __init__(self, residues, centre = [0., 0., 0.], angle = 120, chain = "A"):
        super(VirtualHelix310, self).__init__(residues, centre, angle, chain)


class VirtualHelixPi(VirtualHelixAlpha):

    _MAX_AA_DIST = 1.1
    _RADIUS      = 2.8
    _TYPE        = "I"

    def __init__(self, residues, centre = [0., 0., 0.], angle = 87, chain = "A"):
        super(VirtualHelixPi, self).__init__(residues, centre, angle, chain)


if __name__ == '__main__':
    y = VirtualHelix(20, [5, 5, 5])
    print y.atom_points(1, seq="AAAAAAAAAAAAAAAAAAAAA")
    y.chain = "C"
    y.shift(x=4.9, y=0, z=0.)
    y.shift_to_origin()
    print y.atom_points(1, seq="AAAAAAAAAAAAAAAAAAAAA")

    y = VirtualHelix(20, [0, 0, 0])
    print y.guide_points(1)
    print y.atom_points(1, seq="AAAAAAAAAAAAAAAAAAAAA")
    print y.atom_points(12, seq="AAAAAAAAAAAAAAAAAAAAA")
    y.chain = "B"
    y.tilt_degrees(z_angle = 45)
    print y.guide_points(1)
    print y.atom_points(5, seq="AAAAAAAAAAAAAAAAAAAAA")
    y.chain = "C"
    y.shift(x = 10)
    y.spin_degrees(angle = 100)
    y.invert_direction()
    print y.guide_points(1)
    print y.atom_points(45, seq="AAAAAAAAAAAAAAAAAAAAA")
