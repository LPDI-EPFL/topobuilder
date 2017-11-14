# -*- coding: utf-8 -*-
# @Author: bonet
# @Date:   2016-04-15 18:03:15
# @Last Modified by:   bonet
# @Last Modified time: 2016-05-02 19:35:09

import numpy as np
from VirtualStructure import VirtualStructure as VS


class VirtualBeta(VS):
    # The pleating causes the distance between alpha[i] and alpha[i+2] to be
    # approximately 6, rather than the 7.6 (2 × 3.8) expected from
    # two fully extended trans peptides.
    _ATOMTYPES    = ("N", "CA", "C", "O", "H")
    #_MAX_AA_DIST  = 3.25
    #_ATOM_CA_DIST = {"N": 2.5, "CA": 3.8, "C": 1.43 , "O": 1.55, "H": 2.5}
    #_RADIUS       = {"N": 0.3, "CA": 1.06, "C": 1.082 , "O": 1.082, "H": 0.3}
    #_SHIFT        = {"N": 0.35, "CA": 0., "C": -0.55 , "O": -1.77, "H": 1.35}

    #_ATOM_CA_DIST = {"N": 2.367, "CA": 3.8, "C": 1.405 , "O": 1.655, "H": 2.101}
    #_RADIUS       = {"N": 0.211, "CA": 1.06, "C": -0.305 , "O": -0.171, "H": 0.068}
    #_SHIFT        = {"N": 0.391, "CA": 0., "C": -0.525 , "O": -1.738, "H": 1.354}

    _ATOM_CA_DIST = {"N": 2.6, "CA": 3.8, "C": 1.8 , "O": 1.9, "H": 2.5}
    _RADIUS       = {"N": 0.25, "CA": 1.06, "C": -0.305 , "O": -0.171, "H": 0.068}
    _SHIFT        = {"N": 0.395, "CA": 0., "C": -0.510 , "O": -1.75, "H": 1.35}

    #_ATOM_CA_DIST = {"N": 1.996, "CA": 0., "C": 1.185 , "O": 1.395, "H": 1.772}
    #_RADIUS       = {"N": 0.35, "CA": 1.06, "C": -0.25 , "O": -0.16, "H": 0.08}


    #_RADIUS       = {"N": 0.211, "CA": 1.06, "C": 1.365 , "O": 1.231, "H": 0.211}
    #_RADIUS       = {"N": 1.271, "CA": 1.06, "C": 0.775 , "O": 0.889, "H": 1.128}
    #_RADIUS       = {"N": 0.01, "CA": 0.01, "C": 0.01 , "O": 0.01, "H": 0.01}
    #_SHIFT        = {"N": 0.391, "CA": 0., "C": -0.525 , "O": -1.738, "H": 1.354}

    #_ATOM_CA_DIST = {"N": 2.5, "CA": 0., "C": 1.43 , "O": 1.55, "H": 2.5}
    #_RADIUS       = {"N": 0.032, "CA": 1.06, "C": -0.4 , "O": -0.383, "H": 0.032}
    #_SHIFT        = {"N": 0.2, "CA": 0., "C": -0.5 , "O": -1.65, "H": 1.35}

    # CHOP780202 beta-sheet propensity AAindex (Chou-Fasman, 1978b)
    # TO 0: G -> 0.75; P -> 0.55
    _AA_STAT = [("A", 0.83), ("L", 1.30), ("R", 0.93), ("K", 0.74), ("N", 0.89),
                ("M", 1.05), ("D", 0.54), ("F", 1.38), ("C", 1.19), ("P", 0.00),
                ("Q", 1.10), ("S", 0.75), ("E", 0.37), ("T", 1.19), ("G", 0.00),
                ("W", 1.37), ("H", 0.87), ("Y", 1.47), ("I", 1.60), ("V", 1.70)]

    _TYPE = "E"

    def __init__(self, residues, centre = [0., 0., 0.], chain = "A"):
        super(VirtualBeta, self).__init__(residues, centre, chain)
        self.last_orientation = self._RADIUS["CA"]
        self.atoms = []
        self.atomtypes = []
        for x in range(len(self.points)):
            self.last_orientation *= -1
            for atomtype in self._ATOMTYPES:
                point = np.copy(self.points[x]) + np.array([self._SHIFT[atomtype] * self.last_orientation, self._ATOM_CA_DIST[atomtype], self._RADIUS[atomtype] * self.last_orientation])
                self.atomtypes.append(atomtype)
                self.atoms.append(point)

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
