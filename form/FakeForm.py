# -*- coding: utf-8 -*-
# @Author: bonet
# @Date:   2016-05-01 14:15:09
# @Last Modified by:   bonet
# @Last Modified time: 2016-05-02 23:44:18

import copy


class FakeForm(object):
    """docstring for FakeForm"""
    def __init__(self, sslist):
        self.sslist = copy.deepcopy(sslist)
        self.id     = "_".join(x.desc for x in self.sslist)
        self.turn   = [0 for x in sslist]
        self.do     = 0
        self.edges  = 0
        self.direct = 0
        self.inter  = 0

    def evaluate(self):
        self.edges  = 1 if self._expected_edges() else 0
        self.direct = 1 if self._expected_directions() else 0
        if self.direct == 0:
            self._alternate_regardless()
        self.inter  = 1 if self._expected_intersection() else 0
        evals       = [self.edges, self.direct, self.inter]
        self.do     = sum(evals) == len(evals)

    def to_json(self):
        return { "id": self.id, "do": self.do, "up": self.turn,
                 "obeys": {"edges": self.edges,
                           "directions": self.direct,
                           "intersections": self.inter
                           } }

    def get_ss_by_id(self, name):
        for x in self.sslist:
            if x.desc == name:
                return x

    def _alternate_regardless(self):
        direct = [x.struc.up_is_1() for x in self.sslist]
        directpre = direct[0]
        for x in range(1, len(direct)):
            if direct[x] == directpre:
                self.sslist[x].struc.invert_direction()
                direct[x] = self.sslist[x].struc.up_is_1()
            directpre = direct[x]
        self.turn = direct

    def _expected_directions(self):
        fixed  = [x.static for x in self.sslist]
        direct = [x.struc.up_is_1() for x in self.sslist]
        try:  # start shift at ss 2
            directpre = direct[0]
            fixpre    = fixed[0]
            for x in range(1, len(fixed)):
                if direct[x] == directpre:
                    # 2 consecutive fix in same direction won't work
                    if fixpre and fixed[x]: return False
                    if fixed[x]: raise ValueError
                    self.sslist[x].struc.invert_direction()
                    direct[x] = self.sslist[x].struc.up_is_1()
                directpre = direct[x]
                fixpre    = fixed[x]
        except ValueError:  # start shift at ss 1
            if fixed[0]: return False
            fixpre    = fixed[0]
            self.sslist[0].struc.invert_direction()
            direct[0] = self.sslist[0].struc.up_is_1()
            directpre = direct[0]
            for x in range(1, len(fixed)):
                if direct[x] == directpre:
                    if fixed[x]: return False
                    self.sslist[x].struc.invert_direction()
                    direct[x] = self.sslist[x].struc.up_is_1()
                directpre = direct[x]
                fixpre    = fixed[x]
        self.turn = direct
        return True

    def _expected_edges(self):
        edges = 0
        if self.sslist[0].edge == -1 or self.sslist[-1].edge == -1: return False
        for x in self.sslist:
            if x.edge == 1: edges += 1
        if edges > 2: raise AttributeError("More than two edge structures is not possible")
        if edges != 0 and self.sslist[0].edge + self.sslist[-1].edge != edges: return False
        return True

    def _expected_intersection(self):
        dU = []
        dD = []
        up = 1 if self.sslist[0].struc.goes_up() else -1

        for x in range(len(self.sslist) - 1):
            if up == 1:  dU.append((self.sslist[x], self.sslist[x + 1]))
            if up == -1: dD.append((self.sslist[x], self.sslist[x + 1]))
            up *= -1

        if len(dU) > 1:
            for x in dU:
                for y in dU:
                    if x != y and self._intersection(x, y, "up"):  return False
        if len(dD) > 1:
            for x in dD:
                for y in dD:
                    if x != y and self._intersection(x, y, "down"): return False
        return True

    def _intersection(self, s1, s2, key):
        left   = max(
            min(s1[0].get_x(key), s1[1].get_x(key)),
            min(s2[0].get_x(key), s2[1].get_x(key))
        )
        right  = min(
            max(s1[0].get_x(key), s1[1].get_x(key)),
            max(s2[0].get_x(key), s2[1].get_x(key))
        )
        top    = max(
            min(s1[0].get_z(key), s1[1].get_z(key)),
            min(s2[0].get_z(key), s2[1].get_z(key))
        )
        bottom = min(
            max(s1[0].get_z(key), s1[1].get_z(key)),
            max(s2[0].get_z(key), s2[1].get_z(key))
        )
        if bottom > top and right > left: return True
        return False
