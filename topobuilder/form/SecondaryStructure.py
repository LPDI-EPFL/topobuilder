# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>
.. codeauthor:: Zander Harteveld <zandermilanh@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""

import re
from collections import Iterable
import numpy as np


class SecondaryStructure( object ):
    """docstring for SecondaryStructure"""
    def __init__(self, desc):
        self.desc    = desc["id"]
        self.layer   = ord(re.match('(\w)\d+', self.desc).group(1)) - 64
        self.lcount  = int(re.match('\w(\d+)', self.desc).group(1))
        self.type    = str(desc["type"])
        self.edge    = int(desc["edge"])      if "edge"      in desc else 0
        self.static  = int(desc["static"])    if "static"    in desc else 0
        self.ref     = desc["ref"]            if "ref"       in desc else None
        self._seq    = int(desc["sequence"])  if "sequence"  in desc else None

        self.length  = desc["length"] if "length" in desc else 0
        self.struc   = None

        self.shift_x = desc["shift_x"] if "shift_x" in desc else 0.0
        self.shift_y = desc["shift_y"] if "shift_y" in desc else 0.0
        self.shift_z = desc["shift_z"] if "shift_z" in desc else 0.0

        self.tilt_x = desc["tilt_x"] if "tilt_x" in desc else 0.0
        self.tilt_y = desc["tilt_y"] if "tilt_y" in desc else 0.0
        self.tilt_z = desc["tilt_z"] if "tilt_z" in desc else 0.0

    def get_ref_motif(self):
        return self.ref.split(".")[0]

    def get_ref_segment(self):
        return self.ref.split(".")[1]

    def add_structure(self, struct):
        self.struc = struct

    def get_xyz(self, key):
        e0 = self.struc.edges[0]
        e1 = self.struc.edges[1]
        if key == "up":
            if e0[1] < e1[1]: return e1
            else: return e0
        if key == "down":
            if e0[1] < e1[1]: return e0
            else: return e1

    def get_x(self, key):
        return 0.0 if np.isclose(self.get_xyz(key)[0], 0.0) else self.get_xyz(key)[0]

    def get_y(self, key):
        return 0.0 if np.isclose(self.get_xyz(key)[1], 0.0) else self.get_xyz(key)[1]

    def get_z(self, key):
        return 0.0 if np.isclose(self.get_xyz(key)[2], 0.0) else self.get_xyz(key)[2]

    def twoD_distance( self, other ):
        a = np.array([self.struc.centre[0], self.struc.centre[2]])
        b = np.array([other.struc.centre[0], other.struc.centre[2]])
        return np.linalg.norm(a - b)

    def __eq__( self, other ):
        return self.__hash__() == other.__hash__()

    def __ne__( self, other ):
        return not self.__eq__( other )

    def __lt__( self, other ):
        return self.__hash__() < other.__hash__()

    def __hash__( self ):
        return self.layer * 100 + self.lcount

    def __str__( self ):
        return self.desc

    def __repr__( self ):
        return self.__str__()
