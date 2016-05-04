# -*- coding: utf-8 -*-
# @Author: bonet
# @Date:   2016-04-26 21:35:19
# @Last Modified by:   bonet
# @Last Modified time: 2016-04-26 21:35:27


class Constraint(object):
    """single Constraint"""
    def __init__(self, num1, num2, value, ctype="AtomPair", atm1="CA",
                 atm2="CA", func="HARMONIC", dev=3.0, tag="TAG"):
        self.ctype = ctype
        self.atm1  = atm1
        self.num1  = num1
        self.atm2  = atm2
        self.num2  = num2
        self.func  = func
        self.value = value
        self.dev   = dev
        self.tag   = tag

    def __str__(self):
        atm1 = "{0.atm1} {0.num1}".format(self)
        atm2 = "{0.atm2} {0.num2}".format(self)
        func = "{0.func} {0.value:.2f} {0.dev:.1f} {0.tag}".format(self)
        return "{0.ctype} {1} {2} {3}".format(self, atm1, atm2, func)


class ConstraintSet(object):
    """ConstraintSet"""
    def __init__(self):
        self.constraints = []
        self.constrnsidx = {}

    @staticmethod
    def parse(filename):
        c = ConstraintSet()
        with open(filename) as fd:
            for line in fd:
                l = line.strip().split()
                c.add_constraint(l[2], l[4], l[6], l[0], l[1], l[3], l[5], l[7], l[8])
        return c

    def add_constraint(self, num1, num2, value, ctype="AtomPair", atm1="CA",
                       atm2="CA", func="HARMONIC", dev=3.0, tag="TAG"):
        c = Constraint(int(num1), int(num2), value, ctype, atm1, atm2, func, dev, tag)
        self.constraints.append(c)
        self.constrnsidx.setdefault(int(num1), {})[int(num2)] = c
        self.constrnsidx.setdefault(int(num2), {})[int(num1)] = c

    def has_contact(self, r1, r2):
        return int(r1) in self.constrnsidx and int(r2) in self.constrnsidx[int(r1)]

    def get_contact(self, r1, r2):
        return self.constrnsidx[int(r1)][int(r2)]

    def __getitem__(self, key):
        return self.constraints[key]

    def __len__(self):
        return len(self.constraints)

    def __str__(self):
        text = []
        for x in range(len(self.constraints)):
            text.append(str(self.constraints[x]))
        return "\n".join(text)
