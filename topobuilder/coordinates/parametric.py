# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. class:: ParametricStructure
"""
# Standard Libraries

# External Libraries
import numpy as np
import pandas as pd

# This Library
try:
    from SBI.structure import PDB, ResidueFrame, Frame3D
except ImportError:
    pass


__all__ = ['ParametricStructure']


class ParametricStructure( object ):

    _MONO = None
    _PERIODE = None
    _ROTATION = None

    def __init__( self, indata ):
        """
        """
        if isinstance(indata, Frame3D):
            self.pdb = indata
            self.desc = None
            self.reverse()
        elif isinstance(indata, dict):
            self.pdb = []
            self.desc = indata
            self.build()

    def build( self ):
        """Builds the sturcture, e.g. place the backbone atoms with respect to the
        secondary structure element specified.
        """
        if self._MONO is None or self._PERIODE is None:
            raise NotImplementedError()

        # 1. Locate center point for each residue we need to build
        vector_module = float(self._PERIODE * (self.desc['length'] - 1))
        upper_bound = np.copy(np.array([0., 0., 0.], dtype='float64')) + np.array([0, vector_module / 2, 0])
        points = [np.copy(upper_bound) - np.array([0, self._PERIODE * x, 0]) for x in range(self.desc['length'])]

        # 2. Build. For each point, we build one periode at [0, 0, 0]. Then, we rotate and then shift.
        self.pdb = []
        for i, p in enumerate(points):
            self.pdb.append(ResidueFrame().simple_residue('GLY', self._MONO))
            self.pdb[-1].rotate_degrees(y=self._ROTATION * i)
            self.pdb[-1].translate(p)
            self.pdb[-1].number = i + 1
        self.pdb = pd.concat(self.pdb)

    def reverse( self ):
        """
        """
        raise NotImplementedError()
