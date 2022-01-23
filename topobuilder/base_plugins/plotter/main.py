# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>
.. codeauthor:: Zander Harteveld <zandermilanh@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import List, Union, Optional, Dict
from pathlib import Path
import os

# External Libraries
import matplotlib.pyplot as plt

# This Library
from topobuilder.case import Case
from topobuilder.workflow import Node, NodeOptionsError, NodeDataError
import topobuilder.core as TBcore
from . import plot_types as pts

__all__ = ['plotter']

_PLT_TYPES_ = ['sketchXZ', 'sketchXY']


class plotter( Node ):
    """Create visual representations of the :term:`FORM`.

    .. note::
        Depends on the ``system.image`` configuration option

    .. admonition:: To Developers

        This :class:`.Node` overwrite the execute function to be able to generate
        multi-Case images.

    :param outfile: Selected outfile/outdir for the plot.
    :param prefix: File prefix if outfile is a dir.
    :param plot_types: Select the expected plot (from available options).

    For each plot type, a dict can be provided with extra parameters to control plotting.
    Parameters are defined for :meth:`.plot_case_sketch`

    :raises:
        :NodeOptionsError: On **initialization**. If an unknown plot is required.
        :NodeDataError: On **check**. If the required fields to be executed are not there.
        :NodeDataError: On **execution**. If the requested corrections do not meet the criteria imposed by the
            :term:`FORM`.

    """
    REQUIRED_FIELDS = ('topology.architecture', )
    RETURNED_FIELDS = ()
    PLOT_TYPES = ['sketchXZ', 'sketchXY']
    VERSION = 'v1.0'

    def __init__( self, tag: int,
                  outfile: Optional[Union[str, Path]] = None,
                  prefix: Optional[str] = None,
                  plot_types: Optional[List[str]] = None,
                  **kwargs ):
        super(plotter, self).__init__(tag)

        self.outfile = outfile
        self.prefix = prefix
        self.plot_types = [self.PLOT_TYPES[0], ] if plot_types is None else plot_types
        self.plot_params = {}
        for ptype in self.plot_types:
            self.plot_params.setdefault(ptype, kwargs.pop(ptype, {}))

        if len(set(self.plot_types).difference(self.PLOT_TYPES)) > 0:
            raise NodeOptionsError(f'Unknown plot format. Available are: {",".join(self.PLOT_TYPES)}')

    def single_check( self, dummy: Dict ) -> Dict:
        kase = Case(dummy)

        # Check what it needs
        for itag in self.REQUIRED_FIELDS:
            if kase[itag] is None:
                raise NodeDataError(f'Field "{itag}" is required')

        # Include what keywords it adds (in this instance, nothing)
        return kase.data

    def execute( self, data: List[Dict] ) -> List[Dict]:
        kase = Case(data[0])

        # File management
        if self.outfile is None:
            self.outfile = kase.main_path.joinpath('images').resolve()
            self.outfile.mkdir(parents=True, exist_ok=True)
        if isinstance(self.outfile, str):
            self.outfile = Path(self.outfile).resolve()

        # Get output format
        outformat = TBcore.get_option('system', 'image')

        for ptype in self.plot_types:
            if self.outfile.is_dir():
                self.prefix = self.prefix if self.prefix is not None else ".".join([str(os.getppid()), f'{ptype}'])
                thisoutfile = self.outfile.joinpath(".".join([self.prefix, ptype + outformat]))
            else:
                thisoutfile = Path(str(self.outfile) + '.' + ptype + outformat)
            thisoutfile.parent.mkdir(parents=True, exist_ok=True)
            if not TBcore.get_option('system', 'overwrite') and thisoutfile.is_file():
                self.log.warning(f'Unable to overwrite file {thisoutfile}: Already exists')
                continue

            if not self.plot_params[ptype]:
                fig, ax = getattr(pts, ptype)(self.log, [Case(i) for i in data])
            else:
                fig, ax = getattr(pts, ptype)(self.log, [Case(i) for i in data], self.plot_params[ptype])
            plt.tight_layout()
            plt.savefig(str(thisoutfile), dpi=300)

            self.log.info(f'Creating new image at: {str(thisoutfile)}')
            return data

    def single_execute( self, data: Dict ) -> Dict:
        return data
