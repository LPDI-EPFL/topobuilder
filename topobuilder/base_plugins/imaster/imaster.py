# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>
.. codeauthor:: Zander Harteveld <zandermilanh@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import argparse
from pathlib import Path
from operator import itemgetter
from itertools import cycle

# External Libraries

# This Library
import topobuilder.utils as TButil
from topobuilder.case import Case
from analysis import process_master_geometries
from rstoolbox.io import parse_master_file


def options():
    """
    """
    # Parse Arguments
    parser = argparse.ArgumentParser()

    parser.add_argument('-case', dest='case', action='store', required=True)
    parser.add_argument('-master', dest='master', action='store', required=True)
    parser.add_argument('-present', dest='present', action='store', default=None)
    parser.add_argument('-out', dest='out', action='store', required=True)

    return parser.parse_args()


def main( options ):
    """
    """
    # Load MASTER search data.
    masterdf = parse_master_file(options.master, shift_0=True)

    # Case data
    case = Case(Path(options.case))
    # Get connectivities
    sse = case.connectivities_str[0].split('.')
    # Get flips
    flip = cycle([case['configuration.flip_first'], not case['configuration.flip_first']])
    flip = [next(flip) for _ in range(len(sse))]
    # Select only the present ones.
    present = [sse.index(i) for i in options.present.split('.')]
    sse = list(itemgetter(*present)(sse))
    flip = list(itemgetter(*present)(flip))

    # Geometric properties retrieval
    masterdf = process_master_geometries(masterdf, sse, flip)
    # Output data
    masterdf.to_csv(str(options.out) + '.csv', index=False)


if __name__ == '__main__':
    main(options())
