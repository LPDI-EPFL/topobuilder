# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>
.. codeauthor:: Zander Harteveld <zandermilanh@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import List, Union, Dict
import copy

# External Libraries

# This Library
from topobuilder.case import Case
from topobuilder.workflow import Node, NodeOptionsError, NodeDataError


__all__ = ['nomenclator']


class nomenclator( Node ):
    """Alters the ``configuration.name`` of a :class:`.Case` by adding subnames.

    This affects on the creation of the subfolders where the rest of the :class:`.Pipeline`
    will be executed.

    .. note::
        On **execution**, the plugin will not append new subnames when those already exist. For example,
        if ``configuration.name`` is ``1QYS_experiment1_naive`` and ``subnames=['experiment1', 'naive']``,
        the final ``configuration.name`` will still be ``1QYS_experiment1_naive`` and not
        ``1QYS_experiment1_naive_experiment1_naive``. This is to avoid folder recursion generation when
        re-running a previous :class:`.Pipeline`.

    .. caution::
        There are some keywords that cannot be used as a subname due to them generating their own
        **first level** subfolders. These keywords are ``architecture``, ``connectivity``, ``images``
        and ``summary``. Trying to add one of those terms as subname will generate a :class:`.NodeDataError`
        on **check** time.

    .. admonition:: To Developers

        When developing a new plugin, if it is expected to create new **first level** subfolders, they should
        be listed in the class attribute :attr:`.nomenclator.RESERVED_KEYWORDS`. See more on how to
        :ref:`develop your own plugins <make_plugin>`.

    :param subnames: Subnames that will be added to the :class:`.Case` initial name.

    :raises:
        :NodeOptionsError: On **initialization**. If a reserved key is provided as a subname.
        :NodeDataError: On **check**. If the required fields to be executed are not there.

    """
    RESERVED_KEYWORDS = ['architecture', 'connectivity', 'images', 'summary']
    REQUIRED_FIELDS = ('configuration.name', )
    RETURNED_FIELDS = ()
    VERSION = 'v1.0'

    def __init__( self, tag: int, subnames: Union[List[str], str] ):
        super(nomenclator, self).__init__(tag)

        self.subnames = subnames
        # Unify subnames behaviour for 1 to N
        if not isinstance(self.subnames, list):
            self.subnames = [self.subnames, ]

        # Words used as main folders should not be used for subnaming
        if len(set([x.lower() for x in self.subnames]).intersection(set(self.RESERVED_KEYWORDS))) > 0:
            raise NodeOptionsError(f'Keywords {self.RESERVED_KEYWORDS} cannot be used for subnaming.')

    def single_check( self, dummy: Dict ) -> Dict:
        kase = Case(dummy)

        # Check what it needs
        for itag in self.REQUIRED_FIELDS:
            if kase[itag] is None:
                raise NodeDataError(f'Field "{itag}" is required')

        # Include what keywords it adds (in this instance, nothing)
        return kase.data

    def single_execute( self, data: Dict ) -> Dict:
        kase = Case(data)

        # Check name was not already added.
        sn = copy.deepcopy(self.subnames)
        if kase.name.endswith('_'.join(sn)):
            self.log.notice(f'Seems the subnames {"_".join(sn)} already existed.')
            self.log.notice('Will NOT re-append.')
            return kase

        # Add new names
        oname = kase.name
        sn.insert(0, oname)
        kase.data['configuration']['name'] = '_'.join(sn)

        self.log.debug(f'Renamed case {oname} to {kase.name}')
        return kase.data
