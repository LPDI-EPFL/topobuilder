# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>
.. codeauthor:: Zander Harteveld <zandermilanh@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""

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
