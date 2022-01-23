# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>
.. codeauthor:: Zander Harteveld <zandermilanh@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""

import re

from svgwrite.shapes import Circle, Polygon, Line
from svgwrite.container import Group
from svgwrite.drawing import Drawing
import numpy as np


class Triangle(object):  # Equilater
    def __new__(*args, **kwargs):
        rc = kwargs["rc"]  # circumcircle radius
        del(kwargs["rc"])
        a  = rc / (np.sqrt(3) / 3)
        h  = rc / (np.sqrt(3) / 2)
        center = kwargs["center"]
        del(kwargs["center"])
        rotate = kwargs["rotate"] if "rotate" in kwargs else False
        if "rotate" in kwargs: del(kwargs["rotate"])
        rotate = -1 if rotate else 1
        triangle = Polygon([(0., -rc * rotate),
                            (-a / 2, (h / 2) * rotate),
                            (a / 2, (h / 2) * rotate)], **kwargs)
        triangle.translate(center[0], center[1])
        return triangle


class Cross(object):
    def __new__(*args, **kwargs):
        r = kwargs["r"]
        del(kwargs["r"])
        center = kwargs["center"]
        del(kwargs["center"])
        m = r / 3
        cross = Polygon([(-m, -r), (m, -r), (m, -m), (r, -m), (r, m), (m, m),
                         (m, r), (-m, r), (-m, m), (-r, m), (-r, -m), (-m, -m)],
                         **kwargs)
        cross.translate(center[0], center[1])
        cross.rotate(45, (0., 0.))
        return cross


class VisualForms(object):
    """docstring for VisualForm"""
    def __init__(self, formlist):
        self.forms = formlist

    def make_svg(self, data):
        layers = data["layers"]
        _PROPORTION = 10
        _SHIFT_PROP = _PROPORTION / 1.5
        for f in self.forms:
            structures = [None] * len(f.id.split("_"))
            linkers    = []
            centers    = {}
            order      = {}
            maxW, maxH = 0, 0
            for cx, x in enumerate(f.id.split("_")):
                order[x] = cx
            for cl, l in enumerate(layers):
                for cs, s in enumerate(l):
                    name = f.id + "__" + s["id"]
                    if s["type"] == "H":
                        color  = "cornflowerblue" if "ref" not in s else "gainsboro"
                        shape = Circle( center = (s["shift_x"] * _SHIFT_PROP, s["shift_z"] * _SHIFT_PROP),
                                        r = 2.3 * _PROPORTION, id = name,
                                        fill = color, stroke="black",
                                        stroke_width = "2")
                    elif s["type"] == "E":
                        color  = "indianred" if "ref" not in s else "salmon"
                        rotate = f.get_ss_by_id(s["id"]).struc.goes_down()
                        shape = Triangle(center = (s["shift_x"] * _SHIFT_PROP, s["shift_z"] * _SHIFT_PROP),
                                         rc = 2.3 * _PROPORTION, id = name, fill=color,
                                         stroke = "black", stroke_width = "2",
                                         rotate = rotate)
                    else:
                        color  = "darkgreen" if "ref" not in s else "lightgreen"
                        shape = Cross(center = (s["shift_x"] * _SHIFT_PROP, s["shift_z"] * _SHIFT_PROP), r = 2.3 * _PROPORTION,
                                      fill=color, stroke="black", stroke_width = "2",
                                      id = name)
                    structures[order[s["id"]]] = shape
                    centers[s["id"]] = [s["shift_x"] * _SHIFT_PROP, s["shift_z"] * _SHIFT_PROP]
                    if s["shift_x"] * _SHIFT_PROP > maxW: maxW = s["shift_x"] * _SHIFT_PROP
                    if s["shift_z"] * _SHIFT_PROP > maxH: maxH = s["shift_z"] * _SHIFT_PROP

            for cx, x in enumerate(f.id.split("_")):
                if cx == 0 or cx == len(f.id.split("_")) - 1:
                    init = [0, 0]
                    if (ord(x[0]) - 64) <= len(layers) / 2:
                        init[1] = centers[x][1] - (2.3 * _PROPORTION * 2)
                    else:
                        init[1] = centers[x][1] + (2.3 * _PROPORTION * 2)
                    if int(re.search("(\d+)", x).group(1)) <= len(layers[ord(x[0]) - 65]) / 2:
                        init[0] = centers[x][0] - (2.3 * _PROPORTION * 2)
                    else:
                        init[0] = centers[x][0] + (2.3 * _PROPORTION * 2)
                    if cx == 0:
                        shape = Line(init, centers[x], stroke="darkblue", stroke_width = "4")
                    elif cx == len(f.id.split("_")) - 1:
                        shape = Line(centers[x], init, stroke="darkred", stroke_width = "4")
                    linkers.append(shape)
                for cy, y in enumerate(f.id.split("_")):
                    if cy == cx + 1:
                        shape = Line(centers[x], centers[y], stroke="black", stroke_width = "4")
                        linkers.append(shape)

            # Intercalate
            toplink = []
            dowlink = []
            if f.sslist[0].struc.goes_up():
                dowlink = linkers[0:][::2]
                toplink = linkers[1:][::2]
            else:
                dowlink = linkers[1:][::2]
                toplink = linkers[0:][::2]

            g = Group()
            for x in dowlink:
                g.add(x)
            for x in structures:
                g.add(x)
            for x in toplink:
                g.add(x)
            g.translate(2.3 * _PROPORTION * 3, 2.3 * _PROPORTION * 3)
            d = Drawing(size = (maxW + ((2.3 * _PROPORTION * 3) * 2),
                                maxH + ((2.3 * _PROPORTION * 3) * 2)),
                        id = f.id + "__image")
            d.add(g)
            for x in data["forms"]:
                if x["id"] == f.id:
                    x["svg"] = d.tostring()
