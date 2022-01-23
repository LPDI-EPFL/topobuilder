# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>
.. codeauthor:: Zander Harteveld <zandermilanh@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""

import scipy
import math
import os

from ..RosettaIO.constraints.ConstraintSet import ConstraintSet
from ..RosettaIO.loops.Loops import Loops


class Form(object):
    """docstring for Form"""
    def __init__(self, identifier, sslist, l_linkers):
        self.l_linkers = l_linkers
        self.sslist  = sslist
        self.id      = identifier
        self.seq_str = []
        self.inits   = []
        self.const   = ConstraintSet()
        self.loops   = Loops()
        self.order   = []
        self.insert  = []

    def set_order(self, data):
        order = {}
        for x in range(len(data)):
            order[data[x]] = x + 1
        refs  = filter(None, [x.ref for x in self.sslist])
        self.order = [order[x] for x in refs]

    def make_loops(self):
        for x in range(len(self.sslist)):
            #for i,aa_type in enumerate(self.sslist[x].atomtypes):
                #if aa_type ==
            if self.sslist[x].ref is not None:
                self.loops.add_loop(self.inits[x], self.inits[x] + (len(self.sslist[x].atoms)/4) - 1) # Make cleaner with residue object !
            #if self.sslist[x].ref is not None:
                #self.loops.add_loop(self.inits[x], self.inits[x] + len(self.sslist[x].atoms) - 1)

    def make_constraints(self):
        for x in range(len(self.sslist)):
            y = self.sslist[x]
            p = self.inits[x]
            inner_range = 1 if y.get_type() == 'C' else (2 if y.get_type() == 'E' else 5)
            for r1 in range(1,len(y.atoms),4): # only take Calphas
                for r2 in range(r1 + inner_range*4, len(y.atoms),20): # 1 go to Calphas
                    d = scipy.spatial.distance.euclidean(y.atoms[r1], y.atoms[r2])
                    self.const.add_constraint(num1 = p + r1/4, num2 = p + r2/4, value = d, dev=1.5, tag="INNER")

        for x in range(len(self.sslist)):
            px = self.inits[x]
            sx = self.sslist[x]
            for y in range(x + 1, len(self.sslist)):
                py = self.inits[y]
                sy = self.sslist[y]
                for r1 in range(1,len(sx.atoms),4): # only take Calphas
                    for r2 in range(1,len(sy.atoms),4):
                        d = scipy.spatial.distance.euclidean(sx.atoms[r1], sy.atoms[r2])
                        self.const.add_constraint(num1 = r1/4 + px, num2 = r2/4 + py, value = d, dev=3.0, tag="OUTER")

    # def _check_invert(self):  # TODO: wrong
    #     count1 = 0
    #     for x in range(len(self.sslist)):
    #         if self.sslist[x].ref is not None:
    #             count1 = x
    #     return 0 if count1 % 2 == 1 else 1

    def prepare_coords(self):
        # inv = self._check_invert()
        # for x in range(len(self.sslist)):
        #     if x % 2 == inv:
        #         self.sslist[x].struc.invert_direction()
        if self.l_linkers:
            assert len(self.l_linkers) == len(self.sslist) + 1, "Uppps, did you forget to add the length of the termini (otherwise specify by 'x')."

        if self.l_linkers != None and self.l_linkers[0] != "x":
            if self.l_linkers[0] > 0:
                i = self.l_linkers[0]
                self.inits.append(i)
                for x in range(self.l_linkers[0]):
                    self.seq_str.append(("G", "C", "X"))
            else:# self.l_linkers[0] == 0:
                i = 1
                self.inits.append(i)
        else:
            i = 2
            self.seq_str.append(("G", "C", "X"))
            self.inits.append(i)

        for x in range(len(self.sslist) - 1):
            if self.sslist[x].sequence is None:
                self.sslist[x].create_stat_sequence()
            for xx in self.sslist[x].sequence:
                self.seq_str.append((xx, self.sslist[x].get_type(), "S"))
            i += len(self.sslist[x].sequence)
            #d = scipy.spatial.distance.euclidean(self.sslist[x].atoms[-1], self.sslist[x + 1].atoms[0])
            if self.l_linkers != None and self.l_linkers[x+1] != "x":
                #if (len(self.sslist)-1) < len(self.l_linkers) or self.l_linkers[0] > 0:
                d = self.l_linkers[x+1]
                #else:
                    #d = self.l_linkers[x]
            else:
                d = scipy.spatial.distance.euclidean(self.sslist[x].atoms[-3], self.sslist[x + 1].atoms[1])
                d = int(math.ceil(d / 3.))
            i += d
            for yy in range(d):
                self.seq_str.append(("G", "C", "X"))
            self.inits.append(i)
        if self.sslist[-1].sequence is None:
                self.sslist[-1].create_stat_sequence()
        for xx in self.sslist[-1].sequence:
            self.seq_str.append((xx, self.sslist[-1].get_type(), "S"))
        if self.l_linkers != None and (self.l_linkers[-1] != "x" and self.l_linkers[-1] != 0):
            #if (len(self.sslist)-1) < len(self.l_linkers):
            for x in range(self.l_linkers[-1]):
                self.seq_str.append(("G", "C", "X"))
        else:
            self.seq_str.append(("G", "C", "X"))

        pre=True
        self.insert=[]
        tmp = []
        ins1 = '<Prepend anchor_atom="N" anchor_rsd="{}" connecting_atom="C" jump="false" repeat="1" resname="GLY"/>'
        ins2 = '<Insert anchor_atom="C" anchor_rsd="{}" connecting_atom="N" jump="false" repeat="1" resname="GLY"/>'
        for i, x in enumerate(self.seq_str):
            if x[0] != 'G':
                pre = False
            if x[0] == 'G':
                if not pre:
                    self.insert.append(ins2.format(i))
                else:
                    tmp.append(ins1.format(i + 1))
        self.insert = tmp[::-1] + self.insert


    def to_sequence(self):
        return ">" + self.id + "\n" + "".join([x[0] for x in self.seq_str])

    def to_psipred_ss(self):
        text = []
        text.append("# PSIPRED VFORMAT (PSIPRED V2.6 by David Jones)\n")
        sse = [x[1] for x in self.seq_str]
        seq = [x[0] for x in self.seq_str]
        for i in range(len(sse)):
            pH, pE, pC = 0, 0, 0
            if sse[i] == 'C': pC = 1
            else:
                edge = False
                if self.l_linkers != None and self.l_linkers[-1] == 0:
                    if sse[i] != sse[i - 1] or (sse[i] != sse[i - 2] and sse[i] == 'E'): edge = True
                else:
                    if sse[i] != sse[i - 1] or (sse[i] != sse[i - 2] and sse[i] == 'E'): edge = True
                    if sse[i] != sse[i + 1] or (sse[i] != sse[i + 2] and sse[i] == 'E'): edge = True
                if edge:
                    pC = 0.3
                    if sse[i] == 'E': pE = 0.7
                    else:             pH = 0.7
                else:
                    if sse[i] == 'E': pE = 1
                    else:             pH = 1

            line = "{0:>4} {1} {2}   {3:0.3f}  {4:0.3f}  {5:0.3f}".format(i + 1, seq[i], sse[i], pC, pH, pE)
            text.append(line)
        return '\n'.join(text)

    def to_pdb(self):
        data = []
        # ssdef = []
        for x in range(len(self.sslist)):
            data.append(self.sslist[x].atom_points(atom = self.inits[x]))
        return "\n".join(data)

    def __contains__(self, query):
        for x in self.sslist:
            if x.get_type().upper() == query.upper():
                return True
        return False

    # def to_command(self, ori, chain, targetl, templtl, fastaf, ssefil, constrfl, commfl):

    #     text = []
    #     text.append("-fold_from_loops:target:pose {0}".format(os.path.relpath(ori, os.path.dirname(commfl))))
    #     text.append("-fold_from_loops:target:chain {0}".format(chain))
    #     text.append("-fold_from_loops:target:pdb_count")
    #     text.append("-fold_from_loops:target:loops {0}".format(os.path.relpath(targetl, os.path.dirname(commfl))))
    #     text.append("-fold_from_loops:target:order {0}".format(" ".join([str(x) for x in self.order])))

    #     text.append("-in:file:s {0}".format(os.path.relpath(ori, os.path.dirname(commfl))))
    #     text.append("-in:file:fasta {0}".format(os.path.relpath(fastaf, os.path.dirname(commfl))))
    #     text.append("-in:file:psipred_ss2 {0}".format(os.path.relpath(ssefil, os.path.dirname(commfl))))
    #     text.append("-loops:loop_file {0}".format(os.path.relpath(templtl, os.path.dirname(commfl))))
    #     text.append("-fold_from_loops:scaffold:nopose")
    #     text.append("-fold_from_loops:scaffold:cst_file {0}".format(os.path.relpath(constrfl, os.path.dirname(commfl))))

    #     # text.append("-fold_from_loops:loop_mov_nterm 2")
    #     # text.append("-fold_from_loops:loop_mov_cterm 2")

    #     text.append("-fold_from_loops:hb_srbb 2")

    #     text.append("-fold_from_loops:native_ca_cst")

    #     text.append("-in:ignore_unrecognized_res")
    #     text.append("-run:intermediate_structures")
    #     text.append("-out:nstruct 100")
    #     text.append("-out:file:silent_struct_type binary")
    #     # text.append("-out:prefix {0}".format(ident))
    #     # text.append("-out:file:silent {0}".format(ident))

    #     text.append("-out:overwrite")

    #     return "\n".join(text)
