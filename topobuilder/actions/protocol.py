# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. func:: protocol
"""
# Standard Libraries
from typing import Union, Dict, Optional
from pathlib import Path
import json

# External Libraries
import yaml

# This Library
from topobuilder.case import Case
from topobuilder.workflow import Pipeline, ProtocolIncompatibilityError, EmptyProtocolError

__all__ = ['protocol']


def protocol( case: Union[str, Path, Dict, Case],
              protocol: Optional[Union[str, Path]] = None,
              overwrite: Optional[bool] = False ):
    """
    """
    kase = Case(case)

    protocols = kase['configuration.protocols']
    if protocols is not None:
        if len(protocols) == 1 and not bool(protocols[0]):
            protocols = None
    if protocols is None and protocol is None:
        raise EmptyProtocolError('There are no protocols to run')
    if protocol is not None and protocols is not None:
        raise ProtocolIncompatibilityError('Protocols are provided both through file and in the Case. Pick one.')
    if protocol is not None:
        protocol = str(Path(protocol).resolve())
        try:
            protocols = json.loads("".join([x.strip() for x in open(protocol).readlines()]))
            case_format = 'json'
        except json.JSONDecodeError:
            protocols = yaml.load(open(protocol), Loader=yaml.Loader)
            case_format = 'yaml'

    p = Pipeline(protocols).check(kase.data)
    cases = p.execute([kase.data, ])
    for i, c in enumerate(cases):
        cases[i] = Case(c).assign_protocols(protocols)
        p.log.notice(f'New case file created at {cases[i].write(format=case_format)}')
    return cases
