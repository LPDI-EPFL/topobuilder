# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from pathlib import Path
from typing import Union, List, Dict
import json
import textwrap
import copy

# External Libraries
import yaml

# This Library
import topobuilder.core as TBcore
from topobuilder.case import Case
from topobuilder.utils import IpyExit
import topobuilder

__all__ = ['InteractiveProtocol']


class InteractiveProtocol( object ):
    def __init__( self, case: Union[Case, List[Case]], protocol: Union[Path, Dict] ):
        if not TBcore.get_option('system', 'jupyter'):
            raise topobuilder.interface.InterfaceError('InteractiveProtocol is aimed for IPython only.')

        self.case = None
        self._bckcase = None
        self._load_case(case)

        if isinstance(protocol, Path):
            try:
                protocol = json.loads("".join([x.strip() for x in open(protocol).readlines()]))
            except json.JSONDecodeError:
                protocol = yaml.load(open(protocol), Loader=yaml.Loader)
        self.protocols = protocol
        self.case_count = [0, ] * len(self.protocols)
        self.current = -1

    def reset( self ):
        """Return protocol to the begining and reset :class:`.Case`. Generated data is not deleted.
        """
        self.current = -1
        self.case = copy.deepcopy(self._bckcase)
        self.case_count = [0, ] * len(self.protocols)

    def next( self ):
        """Execute the next plugin on the protocol list.
        """
        next(self)

    def _load_case( self, case ):
        if isinstance(case, Case):
            self.case = [Case(case), ]
            self._bckcase = [Case(case), ]
        elif isinstance(case, list):
            self.case = [Case(x) for x in case]
            self._bckcase = [Case(x) for x in case]

    def __next__( self ):
        if self.current + 1 >= len(self.protocols):
            raise StopIteration('Finished interactive protocols.')
        self.current += 1

        # Get current plugin identifier
        current = copy.deepcopy(self.protocols[self.current])
        name = current.pop('name', None)

        # Execute current plugin
        try:
            self.case = topobuilder.plugin_source.load_plugin(name).apply(self.case, **current, prtid=-1)
            self.case_count[self.current] = len(self.case)
        except IpyExit:  # IpyExit controls questions to the user for interactive behaviour
            self.current -= 1

    def _ipython_display_( self ):
        from IPython.core.display import display, HTML
        from jinja2 import Template

        cols = max([len(x) for x in self.protocols]) * 2
        tplt = Template(textwrap.dedent("""\
            <table class="case_table" style="font-family:monospace;">
                <tr style="border-bottom: 1pt solid black;">
                    <th style="text-align:center;" colspan="{{cols - 1}}"><b>PROTOCOLS</b></th>
                    <th style="text-align:center;"><b>EXECUTED</b></th>
                    <th style="text-align:center;"><b>CASES</b></th>
                <tr>
                {% for item in protocols %}
                <tr>
                    <td style="text-align:left;"><b>{{item.name}}</b></td>
                    {% for key, value in item.items() %}
                        {% if key != 'name' %}
                            <td style="text-align:left;"><b>{{key}}</b></td>
                            <td style="text-align:left;"><i>{{value}}</i></td>
                        {% endif %}
                    {% endfor %}
                    {% for x in range(cols - (item|length * 2)) %}
                        <td></td>
                    {% endfor %}
                   <td style="text-align:center;border-left: 1pt solid black;"><b>
                    {% if loop.index <= current + 1 %}YES{% else %}NO{% endif %}
                   </b></td>
                   <td style="text-align:right;border-left: 1pt solid black;"><b>
                        {{counter[loop.index - 1]}}
                   </b></td>
                </tr>
                {% endfor %}
            </table>""")).render(cols=cols, protocols=self.protocols,
                                 current=self.current, counter=self.case_count)
        display(HTML(tplt))
