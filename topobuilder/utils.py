# -*- coding: utf-8 -*-
# @Author: bonet
# @Date:   2016-04-28 15:03:30
# @Last modified by:   bonet
# @Last modified time: 03-Sep-2019
import os
import shutil
import copy
import numpy as np
from .virtual.VirtualReverse import VirtualReverse
from .virtual.VirtualMaker   import VirtualMaker
from .form.Form import Form
from .RosettaIO.loops.Loops import Loops
from .topoIO import read_all_atom_pdb


def prepare_forms(data, options):
    if data["config"]["status"] >= 4: return

    structures = {}
    for l in data["layers"]:
        for x in l:
            structures[x["id"]] = VirtualMaker(type = x["type"], residues = x["length"])
            structures[x["id"]].name = x["id"]
            if "ref" in x:
                structures[x["id"]].ref = x["ref"]
                group, motif = x["ref"].split(".")
                for y in data["motifs"]:
                    if y["id"] == group:
                        for z in y["segments"]:
                            if z["id"] == motif:
                                structures[x["id"]].add_3AAseq(z["sequence"][::4]) # ONLY TO MAKE IT WORK FOR NOW (4 --> full atom has 4 amino acids)
                                structures[x["id"]].atoms = z["coordinates"]
            else:
                structures[x["id"]].create_stat_sequence()
                structures[x["id"]].tilt_y_degrees(x["tilt_y"])
                structures[x["id"]].tilt_degrees(x["tilt_x"], 0, x["tilt_z"])
            structures[x["id"]].shift(x["shift_x"], x["shift_y"], x["shift_z"])
            structures[x["id"]].remove_movement_memory()

    for x in data["forms"]:
        if x["do"]:
            sslist = []
            refsegs = {}
            for cy, y in enumerate(x["id"].split("_")):
                sslist.append(copy.deepcopy(structures[y]))
                if sslist[-1].up_is_1() != x["up"][cy]:
                    sslist[-1].invert_direction()
                if sslist[-1].ref is not None:
                    refsegs[sslist[-1].ref] = sslist[-1].atoms
            if "l_linkers" not in data["config"]: #data["config"]["l_linkers"] = None
                f = Form(x["id"], sslist, None)
            else:
                f = Form(x["id"], sslist, data["config"]["l_linkers"])
            f.prepare_coords()
            order = []
            for mtf in data["motifs"]:
                for sgm in mtf["segments"]:
                    order.append(mtf["id"] + "." + sgm["id"])
            f.set_order(order)
            f.make_loops()
            f.make_constraints()
            wdir = os.path.join(data["config"]["name"], x["id"])
            if not os.path.isdir(wdir): os.mkdir(wdir)
            seqF, ssF, pdbF, rstF, loopF, cnstF = name_files(wdir)
            with open(seqF, "w") as fd: fd.write(f.to_sequence() + '\n')
            with open(ssF, "w") as fd: fd.write(f.to_psipred_ss())
            with open(pdbF, "w") as fd: fd.write(f.to_pdb())
            with open(loopF, "w") as fd: fd.write(str(f.loops))
            with open(cnstF, "w") as fd: fd.write(str(f.const))

            tmplF, tloF, chain, binders, mEdge = prepare_template(data, wdir, refsegs)


            # CREATE FILES
            tpl = os.path.join(os.path.dirname(os.path.abspath(__file__)), "templates")
            bsh = open(os.path.join(wdir, 'run.sh'), "w")

            # Copy motif file
            shutil.copy(os.path.abspath(data["motifs"][0]["pdbfile"]),
                        os.path.join(wdir, os.path.split(data["motifs"][0]["pdbfile"])[-1]))

            # Copy score file for fragment creation.
            shutil.copy(os.path.join(tpl, 'scores.cfg'), os.path.join(wdir, 'scores.cfg'))
            # Create make_fragments script
            with open(os.path.join(tpl, 'make_fragments.xml')) as fd:
                mkfrags = "".join(fd.readlines())
            mkfrags = mkfrags.format('\n'.join(f.insert))
            with open(os.path.join(wdir, 'make_fragments.xml'), "w") as fd:
                fd.write(mkfrags)
            cmd1 = '{} -parser:protocol make_fragments.xml -in:file:s {} -parser:script_vars vall={}\n'
            bsh.write(cmd1.format(data["config"]["rbin"], os.path.split(pdbF)[-1], data["config"]["vall"]))

            # Create FunFolDes script
            with open(os.path.join(tpl, 'funfoldes.xml')) as fd:
                funfoldes = "".join(fd.readlines())
            funfoldes = funfoldes.format(
                ','.join(['{0[0]}-{0[1]}'.format(_) for _ in f.loops.loops]),
                ','.join(['{0[0]}{1}-{0[1]}{1}'.format(_, mEdge.chain[i_]) for i_, _ in enumerate(mEdge.loops)]),
                mEdge.chain[0],
                os.path.split(data["motifs"][0]["pdbfile"])[-1],
                100
            )
            with open(os.path.join(wdir, 'funfoldes.xml'), "w") as fd:
                fd.write(funfoldes)

            # with open(os.path.join(tpl, "ffl.commands")) as fd:
            #     fflcom = "".join(fd.readlines())
            # fflcom = fflcom.format(tmplF, chain, tloF,
            #                        " ".join([str(j) for j in f.order]),
            #                        os.path.relpath(seqF, wdir),
            #                        os.path.relpath(ssF, wdir),
            #                        os.path.relpath(loopF, wdir),
            #                        os.path.relpath(cnstF, wdir),
            #                        100, 2 if "H" in f else 0, 2 if "E" in f else 0)
            # with open(os.path.join(wdir, "build.commands"), "w") as fd:
            #     fd.write(fflcom)
            #     if len(binders) > 0:
            #         fd.write("-fold_from_loops:target:binders {0}".format(binders))

            with open(os.path.join(tpl, "castor.submiter")) as fd:
                castor = "".join(fd.readlines())
            castor = castor.format(options.user, data["config"]["name"].replace("/", "_") + "_" + x["id"],
                                    200, data["config"]["rbin"])
            with open(os.path.join(wdir, "submiter.sbatch"), "w") as fd:
                fd.write(castor)
            bsh.write('sbatch submiter.sbatch\n')

            bsh.close()

    data["config"]["status"] = 4


def name_files(wdir):
    return (os.path.join(wdir, "sequence.fa"),
            os.path.join(wdir, "structure.ss2"),
            os.path.join(wdir, "sketch.pdb"),
            os.path.join(wdir, "sketch_0001.pdb"),
            os.path.join(wdir, "design.loops"),
            os.path.join(wdir, "constraints.cst"))


def prepare_template(data, wdir, refsegs):
    pdbfile  = ""
    loopfile = ""
    chain    = "A"
    binders  = ""
    if len(data["motifs"]) == 1:
        l = Loops()
        for sgm in data["motifs"][0]["segments"]:
            l.add_loop(sgm["ini"], sgm["end"])
            l.chain.append(data["motifs"][0]["chain"])
        loopfile = os.path.join(data["config"]["name"], "target.loops")
        with open(loopfile, "w") as fd: fd.write(str(l))
        loopfile = os.path.relpath(loopfile, wdir)
        pdbfile = os.path.relpath(data["motifs"][0]["pdbfile"], wdir)
        chain = data["motifs"][0]["chain"]
        binders = data["motifs"][0]["binders"] if "binders" in data["motifs"][0] else ""
    else:
        sections = read_all_atom_pdb(data)
        pdbfile = "target.pdb"
        l = Loops()
        with open(os.path.join(wdir, pdbfile), "w") as fd:
            at, rs = 1, 1
            for s in sections:
                if s.ref is not None:
                    s.align_to(refsegs[s.ref][::4])
                else:
                    s.align_to(refsegs[s.ref])
                fd.write(s.to_pdb(at, rs) + "\n")
                l.add_loop(rs, rs + len(s.guide) - 1)
                at += len(s.atoms)
                rs += len(s.guide)
            fd.write("TER")
            bchain = "B"
            for s in sections:
                s.bind_to_pdb(at, bind_chains = bchain)
                at += sum([len(bx) for bx in s.binat.values()])
                bchain = chr(ord(bchain) + len(s.binat))
        loopfile = "target.loops"
        with open(os.path.join(wdir, loopfile), "w") as fd: fd.write(str(l))

    return pdbfile, loopfile, chain, binders, l


def make_html(data, htmlfile):
    tpl = os.path.join(os.path.dirname(os.path.abspath(__file__)), "templates")
    with open(os.path.join(tpl, "header.html")) as fd:
        header = "".join(fd.readlines())
    with open(os.path.join(tpl, "footer.html")) as fd:
        footer = "".join(fd.readlines())
    with open(os.path.join(tpl, "button-info.html")) as fd:
        button = "".join(fd.readlines())

    formhtml = []
    prestyle = "style=\"border: 1px solid black; text-align: center; padding-bottom: 10px;\""
    for x in data["forms"]:
        conditional = "ng-show=\""
        added       = False
        if not x["obeys"]["edges"]:
            conditional += "!isEdges"
            added = True
        if not x["obeys"]["directions"]:
            if added: conditional += " &&"
            conditional += " !isDir"
            added = True
        if not x["obeys"]["intersections"]:
            if added: conditional += " &&"
            conditional += " !isInt"
            added = True
        if not added: conditional += "true"
        conditional += "\""
        pre  = "<div class=\"col-md-3\" {0} {1}>".format(prestyle, conditional)
        pre += "<div class=\"row\">" + x["id"] + "</div>"
        cnt  = "<div class=\"row\">" + x["svg"] + "</div>"
        info = "<div class=\"row\">"
        info += button.format("success" if x["obeys"]["edges"] else "danger",
                              "success" if x["obeys"]["directions"] else "danger",
                              "success" if x["obeys"]["intersections"] else "danger")
        info += "</div>"
        post = "</div>"
        formhtml.append(pre + cnt + info + post)

    with open(htmlfile, "w") as fd:
        fd.write(header)
        fd.write("\n".join(formhtml))
        fd.write(footer)


def process_motifs(data, options):
    if "motifs" not in data: return
    if data["config"]["status"] >= 2: return

    for mtf in data["motifs"]:
        coord = []
        seq   = []

        # Join Coordinates
        for sgm in mtf["segments"]:
            coord.extend(sgm["coordinates"])
            seq.extend(sgm["sequence"])

        # Reposition
        o = VirtualReverse(mtf["id"], coord, mtf["core"], mtf["lookZ"])
        o.shift_to_origin()
        o.orient_as_origin()

        # Check order, split, fix coords and SS length
        o.fix_direction_to_form(data["layers"], mtf["segments"])
        parts = o.split(mtf["segments"], True, True)
        pdic = {}
        for p in parts:
            pdic[o.name + "." + p.name] = p
        i = -1
        for l in data["layers"]:  # State SS length
            for s in l:
                if "ref" in s and s["ref"] in pdic:
                    s["length"] = len(pdic[s["ref"]].atoms)
                    if i >= 0 and i < len(o.distances):
                        s["shift_x"] = o.distances[i]
                    i += 1
        for c, sgm in enumerate(mtf["segments"]):
            name = o.name + "." + sgm["id"]
            sgm["coordinates"] = map(list, pdic[name].atoms)

    data["config"]["status"] = 2


def form_nomenclator(data):
    letter = "A"
    for cl, l in enumerate(data["layers"]):
        for cs, s in enumerate(l):
            s["id"] = letter + str(cs + 1) + s["type"]
        letter = chr(ord(letter) + 1)
