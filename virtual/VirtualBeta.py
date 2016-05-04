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
    _MAX_AA_DIST  = 3.25
    _RADIUS       = 1.06

    # CHOP780202 beta-sheet propensity AAindex (Chou-Fasman, 1978b)
    # TO 0: G -> 0.75; P -> 0.55
    _AA_STAT = [("A", 0.83), ("L", 1.30), ("R", 0.93), ("K", 0.74), ("N", 0.89),
                ("M", 1.05), ("D", 0.54), ("F", 1.38), ("C", 1.19), ("P", 0.00),
                ("Q", 1.10), ("S", 0.75), ("E", 0.37), ("T", 1.19), ("G", 0.00),
                ("W", 1.37), ("H", 0.87), ("Y", 1.47), ("I", 1.60), ("V", 1.70)]

    _TYPE = "E"

    def __init__(self, residues, centre = [0., 0., 0.], chain = "A"):
        super(VirtualBeta, self).__init__(residues, centre, chain)
        self.last_orientation = self._RADIUS
        for x in range(len(self.points)):
            self.last_orientation *= -1
            point = np.copy(self.points[x]) + np.array([0. , 0., self.last_orientation])
            self.atoms.append(point)


if __name__ == '__main__':
    y = VirtualBeta(12, [0, 0, 0])
    print y.atom_points(1, seq="SSQEALHVTERK")
    y.shift(x=5)
    print y.atom_points(13, seq="SSQEALHVTERK")
    y.shift(x=5)
    y.invert_direction()
    print y.atom_points(25, seq="SSQEALHVTERK")

