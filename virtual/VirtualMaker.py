# -*- coding: utf-8 -*-
# @Author: bonet
# @Date:   2016-04-26 17:25:07
# @Last Modified by:   bonet
# @Last Modified time: 2016-04-28 15:01:03

from .VirtualStructure import VirtualStructure
from .VirtualHelix     import VirtualHelix
from .VirtualBeta      import VirtualBeta


class VirtualMaker(object):
    def __new__(self, *args, **kwargs):
        if "type" in kwargs:
            hid = kwargs["type"].lower()
            if hid in ["alpha", "h", "301", "g", "pi", "i"]: return VirtualHelix(*args, **kwargs)
            del(kwargs["type"])
            if hid == "e": return VirtualBeta(*args, **kwargs)
            else: return VirtualStructure(*args, **kwargs)
