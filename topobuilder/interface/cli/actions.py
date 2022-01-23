# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import argparse
from pathlib import Path

# External Libraries


# This Library
from topobuilder.actions import build, protocol


def cli_build():
    """
    """
    # Parse Arguments
    parser = argparse.ArgumentParser(
        description=cli_build.__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-case', dest='case', action='store',
                        help='Case file.', required=True)
    parser.add_argument('-overwrite', dest='overwrite', action='store_true',
                        default=False, help='Ignore previously existing data. '
                        'Will delete previously existing work folder.')
    options = parser.parse_args()

    build(**vars(options))


def cli_protocol():
    # Parse Arguments
    parser = argparse.ArgumentParser(
        description=cli_build.__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-case', dest='case', action='store',
                        help='Case file.', required=True)
    parser.add_argument('-protocol', dest='protocol', action='store',
                        help='Protocol file.', default=None)
    parser.add_argument('-overwrite', dest='overwrite', action='store_true',
                        default=False, help='Ignore previously existing data. '
                        'Will delete previously existing work folder.')
    options = parser.parse_args()

    options.case = Path(options.case)

    protocol(**vars(options))
