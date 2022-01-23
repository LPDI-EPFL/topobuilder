# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>
.. codeauthor:: Zander Harteveld <zandermilanh@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import List, Dict, Tuple
import itertools
import sys
import textwrap
from string import ascii_uppercase

# External Libraries
import numpy as np

# This Library
from topobuilder.case import Case
import topobuilder.core as TBcore
import topobuilder.utils as TButil

__all__ = ['metadata', 'apply', 'case_apply']


def metadata() -> Dict:
    """Plugin description.

    It includes:

    - ``name``: The plugin identifier.
    - ``Itags``: The metadata tags neccessary to execute.
    - ``Otags``: The metadata tags generated after a successful execution.
    - ``Isngl``: Funtion on the expected input connectivity.
    - ``Osngl``: When :data:`True`, output guarantees single connectivity.
    - ``Ecnd``: If present, apply extra cheching functions.
    """
    def isngl( count: int ) -> bool:
        return count == 0

    def extra( case: Case ) -> Tuple[str, bool]:
        msg = textwrap.dedent("""\
        layer_ranger will delete all previous data from the topology.
        Only empty cases are allowed.""")
        return msg, not case.is_empty

    return {'name': 'nomenclator',
            'Itags': [],
            'Otags': [],
            'Isngl': isngl,
            'Osngl': False,
            'Ecnd': extra}


def apply( cases: List[Case],
           prtid: int,
           ranger: Dict,
           **kwargs ) -> List[Case]:
    """Explore multiple ranges of different secondary structures in each layer.
    """

    TButil.plugin_title(__file__, len(cases))

    new_cases = []
    for i, case in enumerate(cases):
        new_cases.extend(case_apply(case, ranger))

    for i, case in enumerate(new_cases):
        new_cases[i] = case.set_protocol_done(prtid)

    TButil.plugin_case_count(len(new_cases))
    return new_cases


@TButil.plugin_conditions(metadata())
def case_apply( case: Case, ranger: Dict ) -> List[Case]:
    """
    """
    case = Case(case)

    todo = []
    for i, layer in enumerate(ranger):
        sse, rng = layer.split('.')
        rng = [int(x) for x in rng.split('-')]
        if TBcore.get_option('system', 'verbose'):
            sys.stdout.write('Applying {} to {} {} to layer {}\n'.format(*rng, sse, ascii_uppercase[i]))
        todo.append(['{}{}'.format(i, sse) for i in range(rng[0], rng[1] + 1)])

    def empties( lst ):
        return [x for x in lst if not x.startswith('0')]

    new_cases = []
    for architecture in map(empties, itertools.product(*todo)):
        new_cases.append(Case(case.name).add_architecture('.'.join(architecture)))
    return new_cases
