# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import inspect
import re
from pathlib import Path

# External Libraries
import pandas as pd

# This Library
#from topobuilder import plugin_source
import topobuilder

__all__ = ['info_plugins']


def info_plugins() -> pd.DataFrame:
    """
    """
    data = {'name': [], 'description': [], 'argument': [], 'argument_type': [], 'source': []}
    basedir = Path(topobuilder.__file__).parent
    for plugin in plugin_source.list_plugins():
        afunc = plugin_source.load_plugin(plugin).apply
        thisdir = Path(plugin_source.load_plugin(plugin).__file__).parent
        annot = inspect.getfullargspec(afunc)[-1]
        for a in annot:
            if a == 'return':
                continue
            data['name'].append(plugin)
            data['description'].append(afunc.__doc__.split('\n')[0])
            data['argument'].append(a)
            if thisdir == Path().home().joinpath('.topobuilder'):
                data['source'].append('HOMEPLUGIN')
            try:
                if thisdir.relative_to(basedir):
                    data['source'].append('BASEPLUGIN')
            except ValueError:
                data['source'].append('ENVPLUGIN')
            at = str(annot[a]).replace('<class \'', '').replace('\'>', '')
            at = re.sub(r'[a-zA-Z\.]*\.', '', at)
            data['argument_type'].append(at)
    return pd.DataFrame(data).set_index(['source', 'name', 'description', 'argument']).sort_index(level=[0, 1])
