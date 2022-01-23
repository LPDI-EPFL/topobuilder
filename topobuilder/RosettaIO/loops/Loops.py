# -*- coding: utf-8 -*-
# @Author: bonet
# @Date:   2016-04-18 11:10:53
# @Last Modified by:   bonet
# @Last Modified time: 2016-04-18 11:11:35



class Loops(object):
    def __init__(self):
        self.loops = []
        self.chain = []


    def add_loop(self, ini, end):
        self.loops.append((ini, end))

    def __str__(self):
        text = []
        text.append("#LOOP start end cutpoint skip-rate extend")
        for l in self.loops:
            text.append("LOOP {0[0]} {0[1]} 0 0.0 1".format(l))
        return "\n".join(text) + "\n"
