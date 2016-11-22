# -*-
# @project: topobuilder
# @file:    TopobuilderUtils.py
#
# @author: jaume.bonet
# @email:  jaume.bonet@gmail.com
# @url:    jaumebonet.cat
#
# @date:   2016-11-09 16:13:01
#
# @last modified by:   jaume.bonet
# @last modified time: 2016-11-09 17:45:18
#
# -*-
import numpy as np
import sys
import json
from topobuilder.virtual.VirtualReverse import VirtualReverse

def read_motif(pdbfile, in_segments):
    def good_chain(line, chain):
        if not line.startswith("ATOM"): return False
        if line[13:15] != "CA":         return False
        if line[21] not in chain:       return False
        return True

    all_atoms = {}

    for sgm in in_segments["segments"]:
        for i in range(int(sgm["ini"]), int(sgm["end"]) + 1):
            all_atoms.setdefault(str(sgm["chain"]), [])

    with open(pdbfile) as fd:
        for line in fd:
            if good_chain(line, all_atoms.keys()):
                all_atoms[line[21]].append([int(line[22:26]), line[17:20],
                                            np.array([float(line[27:38].strip()),
                                                      float(line[38:46].strip()),
                                                      float(line[46:54].strip())],
                                                     dtype = np.float64)])
    core = {}
    for k in all_atoms:
        for line in all_atoms[k]:
            core.setdefault(k, []).append(line[2])
        core[k] = np.mean(core[k], axis = 0, dtype = np.float64).tolist()

    meancore = []
    for x, sgm in enumerate(in_segments["segments"]):
        for line in all_atoms[sgm["chain"]]:
            if line[0] > int(sgm["end"]):
                in_segments["segments"][x]["core"] = core[sgm["chain"]]
                meancore.append(core[sgm["chain"]])
                break
            if line[0] >= int(sgm["ini"]):
                in_segments["segments"][x].setdefault("sequence", []).append(line[1])
                in_segments["segments"][x].setdefault("coordinates", []).append(line[2])
    in_segments["core"] = np.mean(meancore, axis = 0, dtype = np.float64).tolist()
    print in_segments
    return in_segments

def process_motifs(mtf):

    coord = []
    seq   = []

    # Join Coordinates
    for sgm in mtf["segments"]:
        coord.extend(sgm["coordinates"])
        seq.extend(sgm["sequence"])

    # Reposition
    o = VirtualReverse(mtf["id"], coord, mtf["core"])
    o.shift_to_origin()
    o.orient_as_origin()

    # Check order, split, fix coords and SS length
    # o.fix_direction_to_form(data["layers"], mtf["segments"])
    parts = o.split(mtf["segments"], False, False)
    print o.atom_points()
    for cp, p in enumerate(parts):
        mtf["segments"][cp]["center"] = p.centre.tolist()
        del(mtf["segments"][cp]["core"])
    for m in mtf["segments"]:
        for ci, i in enumerate(m["coordinates"]):
            m["coordinates"][ci] = m["coordinates"][ci].tolist()
    del(mtf["core"])
    return mtf

def byteify(input):
    if isinstance(input, dict):
        return {byteify(key): byteify(value)
                for key, value in input.iteritems()}
    elif isinstance(input, list):
        return [byteify(element) for element in input]
    elif isinstance(input, unicode):
        return input.encode('utf-8')
    else:
        return input

# data = '{"segments": [{"end": 69, "chain": "F", "split": "", "$$hashKey": "object:66", "type": "X", "id": "loop", "ini": 62}, {"end": 212, "chain": "F", "split": "", "$$hashKey": "object:76", "type": "H", "id": "helix", "ini": 196}], "id": "d25"}'
# file = '/Users/bonet/SandBox/SketchTopoBenchInteractive/test/4jhw.pdb'

# print json.dumps(process_motifs(read_motif(file, byteify(json.loads(data)))))