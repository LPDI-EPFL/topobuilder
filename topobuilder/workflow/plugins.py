# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import os
from pathlib import Path
from argparse import Namespace

# External Libraries
from pluginbase import PluginBase

# This Library

__all__ = ['TBplugins']


def find_all_plugins():
    """List all possible locations of *TopoBuilder* Plugins.

    By default, this includes 2 main locations:
        * ``base_plugins`` is found inside the *TopoBuilder* distribution.
        * ``$TOPOBULDERPLUGINS`` environmental variable can be used to specify a
            development folder for plugins.
    """

    tp_plugin_dirs = [str(Path(__file__).parent.parent.joinpath('base_plugins')), ]
    try:
        tp_plugin_dirs.append(os.environ['TOPOBULDERPLUGINS'])
    except KeyError:
        pass
    return tp_plugin_dirs


plugin_base = PluginBase(package='topobuilder.plugins')
TBplugins = Namespace(base=plugin_base, source=plugin_base.make_plugin_source(searchpath=find_all_plugins()))
