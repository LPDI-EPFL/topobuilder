# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>
.. codeauthor:: Zander Harteveld <zandermilanh@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""

import copy
import numpy as np
import scipy.spatial
from VirtualStructure import VirtualStructure as VS
from transforms3d.euler import euler2mat

DEFAULT_MODULE = 20


class VirtualReverse(VS):
    """docstring for VirtualReverse"""
    def __init__(self, name, coordinates, coreguide, lookZ=1, chain = "A"):
        coordinates = np.array(coordinates, dtype = np.float64)
        residues = len(coordinates)
        centre = np.mean(coordinates, axis = 0, dtype = np.float64).tolist()
        super(VirtualReverse, self).__init__(residues, centre, chain)
        self.name  = name
        self.atoms = coordinates
        self.get_center_guide_point(coreguide)

        self.MajAx = None  # Point for Major Axis: Line(self.centre, self.MajAx)
        self.MinAx = None  # Point for Minor Axis: Line(self.centre, self.MinAx)
        self.CrsAx = None  # Point for Cross Axis: Line(self.centre, self.CrsAx)
        self.EigVc = None  # EigenVector of the object

        self.distances = []

        self.axis(lookZ)

    # EXTRA SHIFTS
    def shift(self, x = 0, y = 0, z = 0):
        t = np.array(x) if len(x) == 3 else np.array([x, y, z])
        super(VirtualReverse, self).shift(t)
        self.MajAx += t
        self.MinAx += t
        self.CrsAx += t

    # EXTRA ROTATIONS
    def apply_matrix(self, R):
        super(VirtualReverse, self).apply_matrix(R)
        self.MajAx = np.dot(self.MajAx, R)
        self.MinAx = np.dot(self.MinAx, R)
        self.CrsAx = np.dot(self.CrsAx, R)

    # CALCULATE AXIS
    def axis(self, lookZ):
        return self.major_axis(), self.minor_axis(), self.cross_axis(lookZ)

    def major_axis(self):
        if self.MajAx is None:
            self.MajAx = self._calculate_axis(2, 0)
        return self.MajAx

    def minor_axis(self):
        if self.MinAx is None:
            self.MinAx = self._calculate_axis(1, 0)
        return self.MinAx

    def cross_axis(self, lookZ):
        if self.CrsAx is None:
            self.CrsAx = self._calculate_axis(0, 0)

            slope = self.CrsAx - self.centre
            slope = slope / np.linalg.norm(slope)
            icross = self.centre + (-DEFAULT_MODULE) * slope
            a = scipy.spatial.distance.euclidean(self.CrsAx, self.points[0])
            b = scipy.spatial.distance.euclidean(icross, self.points[0])
            if b > a and lookZ == 1: self.CrsAx = icross

        return self.CrsAx

    def _calculate_axis(self, axis, module = 0):
        if not module: module = DEFAULT_MODULE
        if self.EigVc is None:
            A  = np.asmatrix(np.zeros((3, 3)))
            P  = self.atoms - self.centre
            for p in P:
                r = np.asmatrix(p)
                A += r.transpose() * r
            val, self.EigVc = np.linalg.eigh(A)
        t = np.asarray(self.EigVc[:, axis]).reshape(3)
        return self.centre + module * t

    # REORIENT
    def orient_as_origin(self):
        fixpos = self.shift_to_origin()
        x_angle, y_angle, z_angle = self._angles_to_origin()
        self.tilt_radiants(z_angle, y_angle, x_angle)
        self.shift(fixpos)
        self.remove_movement_memory()

    def _angles_to_origin(self):
        assert np.allclose(self.centre, np.array([0, 0, 0]))

        yxz = np.copy(self.axis(None))
        yy  = [0., 1., 0.]
        zz  = [0., 0., 1.]
        # AXIS Y: OVER PLANE Y, X
        p = [yxz[0][0], yxz[0][1], 0.] / np.linalg.norm([yxz[0][0], yxz[0][1], 0.])
        b = np.arccos(np.dot(p, yy)) * -(yxz[0][0] / np.abs(yxz[0][0]))
        Rx = euler2mat(0, 0, b, "sxyz")
        yxz = np.dot(yxz, Rx)
        # AXIS Y: OVER PLANE Y, Z
        p = [0., yxz[0][1], yxz[0][2]] / np.linalg.norm([0., yxz[0][1], yxz[0][2]])
        a = np.arccos(np.dot(p, yy)) * (yxz[0][2] / np.abs(yxz[0][2]))
        Rz = euler2mat(a, 0, 0, "sxyz")
        yxz = np.dot(yxz, Rz)
        # AXIS Z: OVER PLANE Z, X
        p = [yxz[2][0], 0., yxz[2][2]] / np.linalg.norm([yxz[2][0], 0., yxz[2][2]])
        c = np.arccos(np.dot(p, zz)) * (yxz[2][0] / np.abs(yxz[2][0]))
        return b, c, a

    # EXTRA FUNCTIONS
    def get_center_guide_point(self, point):
        self.points[0] = np.array(point, dtype = np.float64)

    def fix_direction_to_form(self, form, segments):
        if len(segments) == 1: return

        parts = self.split(segments, True)
        thisorder = []  # left to right
        for cx, x in enumerate(parts):
            n = str(self.name + "." + x.name)
            if cx == 0: thisorder.append(n)
            else:
                if x.centre[0] < parts[cx - 1].centre[0]: thisorder.insert(-1, n)
                else: thisorder.append(n)

        expected = []  # left to right
        for l in form:
            for ss in l:
                if "ref" in ss and ss["ref"] in thisorder:
                    ss["static"] = 1
                    expected.append(str(ss["ref"]))

        if thisorder == expected:  # Same order
            pass  # Nothing to do
        elif list(reversed(thisorder)) == expected:  # Inverted order
            self.tilt_z_degrees(180)
        else:
            raise Exception("Unable to fit motif in form")

    # SPLIT
    def split(self, segments, correlate_center = False, translate_center = False):
        # if len(segments) == 1:
        #     return copy.deepcopy(self)
        newobj = []
        i = 0
        for sgm in segments:
            lenseg = len(sgm["coordinates"])
            coords = self.atoms[i:i + lenseg]
            newobj.append(VirtualReverse(sgm["id"], coords, self.points[0], 0))
            if correlate_center:
                diff = newobj[-1].centre - self.centre
                diff[0] = 0
                newobj[-1].centre -= diff
            i += lenseg
        for n1 in range(len(newobj)):
            for n2 in range(n1 + 1, len(newobj)):
                self.distances.append(newobj[n1].center_distance_to(newobj[n2]))
        if translate_center:
            for n in newobj:
                n.shift_to_origin()
        return newobj

    # PRINT
    def guide_points(self, atom = 1):
        data = []
        data.append(self._STRING_X.format(atom + 0, self.centre, self.chain))
        data.append(self._STRING_X.format(atom + 1, self.MajAx, self.chain))
        data.append(self._STRING_X.format(atom + 2, self.MinAx, self.chain))
        data.append(self._STRING_X.format(atom + 3, self.CrsAx, self.chain))
        return "\n".join(data)
