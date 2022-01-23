# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import os
import sys
from typing import Union, List, Dict
from functools import wraps

# External Libraries
import colorama as cl

# This Library
import topobuilder.core as TBcore

cl.init(autoreset=True)


def plugin_title( plugin_path: str, cases: int ):
    """Print on-screen the plugin's name.

    :param str plugin_path: ``__file__`` of the plugin's main file.
    :param int cases: Number of cases to which the plugin will be applied.
    """
    name = os.path.basename(os.path.dirname(plugin_path)).upper()
    name = ' '.join(['  TB PLUGIN:', name, ' '])
    bord = ''.join(['_', ] * len(name))

    sys.stdout.write('\n')
    sys.stdout.write(cl.Style.BRIGHT + bord + '\n')
    sys.stdout.write(cl.Style.BRIGHT + name + '\n')
    sys.stdout.write(cl.Style.BRIGHT + bord + '\n')

    plugin_case_count(cases, 'i')


def plugin_case_count( cases: int, io='o' ):
    """Print on-screen the current number of cases.

    :param int cases: Current number of cases.
    :param str io: If ``i``: input. If ``o``: output
    """
    if TBcore.get_option('system', 'verbose'):
        io = 'batch applied to' if io == 'i' else 'generated a total of'
        sys.stdout.write(cl.Style.DIM + '* {} {:03d} cases\n\n'.format(io, cases))


def plugin_warning( text: str ):
    """Format a warning and manage follow up behaviour.

    :param str text: Warning text.
    """
    sys.stdout.write(cl.Fore.GREEN + str(text) + '\n')
    if TBcore.get_option('system', 'strict'):
        sys.exit(-1)


def plugin_filemaker( text: str ):
    """Highlight the creation of new files.

    :param str text: File creation text.
    """
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write(cl.Fore.WHITE + cl.Back.BLUE + 'NEW FILE: ' + str(text))
        sys.stdout.write(cl.Style.RESET_ALL + '\n')


def plugin_imagemaker( text: str ):
    """Highlight the creation of new plot files.

    :param str text: File creation text.
    """
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write(cl.Fore.GREEN + cl.Back.BLUE + 'NEW PLOT: ' + str(text))
        sys.stdout.write(cl.Style.RESET_ALL + '\n')


def plugin_filereader( text: str ):
    """Highlight the loading of files.

    :param str text: File reading text.
    """
    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write(cl.Fore.BLUE + cl.Back.WHITE + 'READ FILE: ' + str(text))
        sys.stdout.write(cl.Style.RESET_ALL + '\n')


def plugin_bash( text: Union[List[List], List[str]] ):
    """Show prepared bash statements to execute.

    :param text: List of statements to show.
    """
    if TBcore.get_option('system', 'verbose'):
        if isinstance(text[0], str):
            sys.stdout.write(cl.Fore.RED + cl.Back.CYAN + 'BASH: ' + ' '.join([str(x) for x in text]))
        else:
            for cmd in text:
                sys.stdout.write(cl.Fore.RED + cl.Back.CYAN + 'BASH: ' + ' '.join([str(x) for x in cmd]) + '\n')
        sys.stdout.write(cl.Style.RESET_ALL + '\n')


def plugin_conditions( metadata: Dict ):
    """Check if the :class:`.Case` has reach the appropiate data content to be executed.

    :param dict metadata: Plugin requirements defined by its ``metadata`` function.
    """
    def wrap(func):
        @wraps(func)
        def wrapper( *args, **kwargs ):
            # Retrieve the Case.
            case = kwargs.get('case', args[0])

            if TBcore.get_option('system', 'debug'):
                sys.stdout.write(cl.Style.DIM + 'Checking viability of plugin {}\n'.format(metadata['name']))

            # Check connectivities
            if not metadata['Isngl'](case.connectivity_count):
                raise PluginOrderError('Plugin {} can only be applied to one connectivity.'.format(metadata['name']))

            # Check metadata
            for tag in metadata['Itags']:
                if case['metadata.{}'.format(tag)] is None:
                    raise PluginOrderError('Missing info from metadata.{}'.format(tag))

            # Extras
            if 'Ecnd' in metadata:
                msg, err = metadata['Ecnd'](case)
                if err:
                    raise PluginOrderError(msg)

            return func(*args, **kwargs)
        return wrapper
    return wrap


class PluginOrderError( Exception ):
    """Raised when plugins expect data from other plugins and do not find it.
    """
