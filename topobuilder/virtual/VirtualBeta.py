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


class VirtualBeta(VS):
    # The pleating causes the distance between alpha[i] and alpha[i+2] to be
    # approximately 6, rather than the 7.6 (2 × 3.8) expected from
    # two fully extended trans peptides.
    _ATOMTYPES    = ("N", "CA", "C", "O")#, "H")

    _ATOM_CA_DIST = {"N": 5.000, "CA": 3.800, "C": 2.600, "O": 2.500} #"H": 2.5}
    _RADIUS       = {"N": 0.300, "CA": 1.100, "C": 0.300, "O": 0.200} #"H": 0.068}
    _SHIFT        = {"N": 0.400, "CA": 0.000, "C": 0.500, "O": 1.900} #"H": 1.35}

    # CHOP780202 beta-sheet propensity AAindex (Chou-Fasman, 1978b)
    # TO 0: G -> 0.75; P -> 0.55; C -> 1.19
    _AA_STAT = [("A", 0.83), ("L", 1.30), ("R", 0.93), ("K", 0.74), ("N", 0.89),
                ("M", 1.05), ("D", 0.54), ("F", 1.38), ("C", 0.00), ("P", 0.00),
                ("Q", 1.10), ("S", 0.75), ("E", 0.37), ("T", 1.19), ("G", 0.00),
                ("W", 1.37), ("H", 0.87), ("Y", 1.47), ("I", 1.60), ("V", 1.70)]

    _TYPE = "E"

    def __init__(self, residues, centre = [0., 0., 0.], chain = "A"):
        super(VirtualBeta, self).__init__(residues, centre, chain)
        self.last_orientation = self._RADIUS["CA"]
        self.atoms = []
        self.atomtypes = []
        self.residuenumbers = []
        self.ca_atoms = []
        count = 0
        for x in range(len(self.points)):
            count += 1
            self.last_orientation *= -1
            #ca_point = np.copy(self.points[x]) + np.array([0. , 0., self.last_orientation])
            #self.atoms.append(ca_point)
            for atomtype in self._ATOMTYPES:
                points = np.copy(self.points[x]) + np.array([self._SHIFT[atomtype] * self.last_orientation, self._ATOM_CA_DIST[atomtype] - self._ATOM_CA_DIST["CA"], self._RADIUS[atomtype] * self.last_orientation])
                self.atomtypes.append(atomtype)
                self.atoms.append(points)
                #if (1 + x)%len(self._ATOMTYPE)==0:
                    #count += 1
                #self.residuenumbers.append(1 + x + count)
                self.residuenumbers.append(count)
                if atomtype == "CA":
                    self.ca_atoms.append(points)

if __name__ == '__main__':
    y = VirtualBeta(16, [0., 0., 0.])

    y.shift_to_origin()
    print y.atom_points(2, seq="AAAAAAAAAAAAAAAA")

    y.shift_to_origin()
    y.shift(x=4.9, y=0, z=0.)
    print y.atom_points(17, seq="AAAAAAAAAAAAAAAA")

    y.shift_to_origin()
    y.shift(x=4.9, y=1., z=0.)
    y.invert_direction()
    print y.atom_points(31, seq="AAAAAAAAAAAAAAAAA")
