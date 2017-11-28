# -*- coding: utf-8 -*-
# @Author: bonet
# @Date:   2016-05-01 12:31:37
# @Last Modified by:   bonet
# @Last Modified time: 2016-05-02 16:58:11
import networkx as nx
import numpy as np
import copy
import sys
import os

from ..virtual.VirtualMaker import VirtualMaker

from .FakeForm import FakeForm
from .Form import Form
from .VisualForms import VisualForms
from .SecondaryStructure import SecondaryStructure as SS

class FormFabric(object):
    """docstring for FormFabric"""
    def __init__(self):
        pass

    def build(self, data, options):
        if data["config"]["status"] >= 3: return

        _DEF_Z_DISTANCE = data["config"]["default_z"] if "default_z" in data["config"] else 11.
        _DEF_X_DISTANCE = {"H": data["config"]["default_x_h"] if "default_x_h" in data["config"] else 11.,
                           "E": data["config"]["default_x_e"] if "default_x_e" in data["config"] else 5.}
        #_DEF_Z_TILT = data["config"]["default_z_tilt"] if "default_z_tilt" in data["config"] else 0.
        _LINK_DISTANCE  = (np.sqrt(2 * (_DEF_Z_DISTANCE * _DEF_Z_DISTANCE))) + 2.0
        _LINK_DISTANCE  = data["config"]["link_dist"] if "link_dist" in data["config"] else _LINK_DISTANCE

        layers    = []
        shapelsit = []
        for x in range(len(data["layers"])):
            layers.append([])
            width = 0
            for y in range(len(data["layers"][x])):
                ss   = data["layers"][x][y]
                if "shift_y" not in ss: ss["shift_y"] = 0.0
                if "shift_z" not in ss: ss["shift_z"] = _DEF_Z_DISTANCE * x
                else: ss["shift_z"] += (_DEF_Z_DISTANCE * x)
                if "shift_x" not in ss:
                    xdist = _DEF_X_DISTANCE["E"] if ss["type"] == "E" else _DEF_X_DISTANCE["H"]
                    xdist = 0.0 if y == 0 else xdist
                    ss["shift_x"] = width + xdist
                else:
                    ss["shift_x"] += width
                if "tilt_x" not in ss: ss["tilt_x"] = 0.0
                if "tilt_y" not in ss: ss["tilt_y"] = 0.0
                if "tilt_z" not in ss: ss["tilt_z"] = 0.0
                #if "tilt_z_layer" not in ss:
                    #ss["tilt_z"] = ztilt_layer
                #else:
                    #ztilt_layer = ss["tilt_z_layer"]
                secstr = SS(ss)
                vs = VirtualMaker(ss["length"], [0., 0., 0.], type=ss["type"])
                if secstr.ref is not None:
                    ref = None
                    for mtf in data["motifs"]:
                        if mtf["id"] == secstr.get_ref_motif():
                            for z in mtf["segments"]:
                                if z["id"] == secstr.get_ref_segment():
                                    ref = z
                                    break
                    vs.atoms = ref["coordinates"]
                else:
                    vs.tilt_y_degrees(ss["tilt_y"])
                    vs.tilt_degrees(ss["tilt_x"], 0, ss["tilt_z"])
                xdist = _DEF_X_DISTANCE["E"] if vs.get_type() == "E" else _DEF_X_DISTANCE["H"]
                vs.shift(ss["shift_x"], ss["shift_y"], ss["shift_z"])
                width = vs.centre[0]
                secstr.add_structure(vs)
                layers[-1].append(secstr)
                shapelsit.append(vs)

        shapeForm = Form("shapesketch", shapelsit)
        shapeForm.prepare_coords()
        with open(os.path.join(options.outdir, "shapesketch.pdb"), "w") as fd:
            fd.write(shapeForm.to_pdb())
        if options.shape: sys.exit(1)

        print "\tevaluating possible combinations"
        forms = _create_forms(layers, _LINK_DISTANCE)
        print "\tforms created:", str(len(forms))

        # EVALUATE AND SAVE FORMS FOR CHECKPOINT
        data.setdefault("forms", [])
        okforms = []
        for _, f in enumerate(forms):
            f.evaluate()
            if f.do: okforms.append(f)
            if f.do or not options.hurry:
                data["forms"].append(f.to_json())
            if _ > 0 and _ % 100 == 0:
                print "\t\t{0} out of {1} evaluated ({2} ok)".format(_, len(forms), len(okforms))
        print "\t\t{0} evaluated ({1} ok)".format(len(forms), len(okforms))

        # GRAPHIC REPRESENTATIONS
        vs = VisualForms(okforms if options.hurry else forms)
        vs.make_svg(data)

        data["config"]["status"] = 3

def _create_forms( layers, distance ):
    G = _create_graph(layers, distance)
    # path_length = len( G.nodes() ) - 1
    forms = []

    for _1, n1 in enumerate(G.nodes()):
        for _2, n2 in enumerate(G.nodes()):
            if _1 < _2:
                forms.extend(_search_paths(G, n1, n2))
                # print "\t\t", n1.desc, "-->", n2.desc
                # for path in nx.all_simple_paths(G, n1, n2):
                #     if len(path) == nx.number_of_nodes(G):
                #         f = FakeForm(copy.deepcopy(path))
                #         forms.append(f)
                #         path.reverse()
                #         f = FakeForm(copy.deepcopy(path))
                #         forms.append(f)

    # for node in G.nodes():
    #     print node.desc
    #     for path in _find_paths(G, node, path_length):
    #         f = FakeForm(copy.deepcopy(path))
    #         forms.append(f)
    return forms

def _search_paths(G, n1, n2):
    forms = []
    print "\t\t", n1.desc, "-->", n2.desc
    for path in nx.all_simple_paths(G, n1, n2):
        if len(path) == nx.number_of_nodes(G):
            f = FakeForm(copy.deepcopy(path))
            forms.append(f)
            path.reverse()
            f = FakeForm(copy.deepcopy(path))
            forms.append(f)
    print "\t\t\t", len(forms), "folds obtained"
    return forms

def _create_graph( layers, distance ):
    G = nx.Graph()
    for lyr1 in range(len(layers)):
        for lyr2 in range(len(layers)):
            if abs(lyr1 - lyr2) <= 1:  # Only consecutive layers
                for col1 in range(len(layers[lyr1])):
                    for col2 in range(len(layers[lyr2])):
                        if abs(col1 - col2) <= 1:  # Only consecutive columns
                            G.add_edge(layers[lyr1][col1],
                                       layers[lyr2][col2], object = SS)
    for lyr1 in layers:
        for lyr2 in layers:
            for sse1 in lyr1:
                for sse2 in lyr2:
                    if sse1 != sse2 and sse1.twoD_distance(sse2) <= distance:
                        G.add_edge(sse1, sse2)

    return G

def _find_paths( G, u, n ):
    if n == 0: return [[u]]
    paths = [[u] + path for neighbor in G.neighbors(u) for path in _find_paths(G, neighbor, n - 1) if u not in path]
    return paths
