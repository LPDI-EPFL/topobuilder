# -*- coding: utf-8 -*-
# @Author: bonet
# @Date:   2016-05-01 12:31:37
# @Last Modified by:   bonet
# @Last Modified time: 2016-05-02 16:58:11
import networkx as nx
import copy

from ..virtual.VirtualMaker import VirtualMaker

from .FakeForm import FakeForm
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

        layers = []
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
                    vs.tilt_degrees(ss["tilt_x"], ss["tilt_y"], ss["tilt_z"])
                xdist = _DEF_X_DISTANCE["E"] if vs.get_type() == "E" else _DEF_X_DISTANCE["H"]
                vs.shift(ss["shift_x"], ss["shift_y"], ss["shift_z"])
                width = vs.centre[0]
                secstr.add_structure(vs)
                layers[-1].append(secstr)

        forms = self._create_forms(layers)

        # EVALUATE AND SAVE FORMS FOR CHECKPOINT
        data.setdefault("forms", [])
        for f in forms:
            f.evaluate()
            data["forms"].append(f.to_json())

        # GRAPHIC REPRESENTATIONS
        vs = VisualForms(forms)
        vs.make_svg(data)

        data["config"]["status"] = 3

    def _create_forms( self, layers ):
        G = self._create_graph(layers)
        path_length = len( G.nodes() ) - 1
        forms = []
        for node in G.nodes():
            for path in self._find_paths(G, node, path_length):
                f = FakeForm(copy.deepcopy(path))
                forms.append(f)
        return forms

    def _create_graph( self, layers ):
        G = nx.Graph()
        for lyr1 in range(len(layers)):
            for lyr2 in range(len(layers)):
                if abs(lyr1 - lyr2) <= 1:  # Only consecutive layers
                    for col1 in range(len(layers[lyr1])):
                        for col2 in range(len(layers[lyr2])):
                            if abs(col1 - col2) <= 1:  # Only consecutive columns
                                G.add_edge(layers[lyr1][col1],
                                           layers[lyr2][col2], object = SS)
        # for x in layers:
        #     for sse1 in x:
        #         for sse2 in x:
        #             if sse1 < sse2:
        #                 G.add_edge( sse1 , sse2, object=SS )
        # for lyr1 in range(len(layers)):
        #     for lyr2 in range(len(layers)):
        #         if abs(lyr1 - lyr2) == 1:  # Only consecutive layers
        #             for sse1 in layers[lyr1]:
        #                 for sse2 in layers[lyr2]:
        #                     G.add_edge( sse1 , sse2, object=SS )
        return G

    def _find_paths( self, G, u, n ):
        if n == 0: return [[u]]
        paths = [[u] + path for neighbor in G.neighbors(u) for path in self._find_paths(G, neighbor, n - 1) if u not in path]
        return paths
