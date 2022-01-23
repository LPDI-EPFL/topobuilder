# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>
.. codeauthor:: Zander Harteveld <zandermilanh@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
from typing import List, Tuple
import sys
import copy

# External Libraries
from logbook import Logger
import matplotlib.pyplot as plt

# This Library
import topobuilder.core as TBcore
from topobuilder.case import Case, plot_case_sketch, plot_case_sketch_vertical

__all__ = ['sketchXZ', 'sketchXY']


def sketchXZ( log: Logger, cases: List[Case], **kwargs ) -> Tuple[plt.Figure, List[plt.Axes]]:
    """
    """
    grid = _calculate_grid(cases, **kwargs)

    if TBcore.get_option('system', 'verbose'):
        sys.stdout.write('Generating an image grid of: {0}x{1}\n'.format(grid[0], grid[1]))

    fsize = (kwargs.pop('width', 7.5 * grid[1]),
             kwargs.pop('hight', 7.5 * grid[0]))
    fig = plt.figure(figsize=fsize)
    axs = []
    ylim, xlim = [0, 0], [0, 0]
    for i, case in enumerate(cases):
        position = (int(i / grid[1]), i % grid[1])
        title = '{0}_{1:03d}'.format(case['configuration.name'], i + 1)
        log.info(f'Showing {title}-{case.architecture_str} in position: {position[0]}x{position[1]}\n')

        ax = plt.subplot2grid(grid, position, fig=fig)
        axs.append(ax)
        plot_case_sketch(case, ax,
                         kwargs.pop('connections', False),
                         kwargs.pop('beta_fill', 'red'),
                         kwargs.pop('beta_edge', 'black'),
                         kwargs.pop('alpha_fill', 'blue'),
                         kwargs.pop('alpha_edge', 'black'),
                         kwargs.pop('connection_edge', None))
        ax.set_title(title)
        cy = ax.get_ylim()
        cx = ax.get_xlim()
        ylim = [ylim[0] if cy[0] < ylim[0] else cy[0],
                ylim[1] if cy[1] > ylim[1] else cy[1]]
        xlim = [xlim[0] if cx[0] > xlim[0] else cx[0],
                xlim[1] if cx[1] < xlim[1] else cx[1]]

    for ax in axs:
        ax.set_ylim(ylim[0], ylim[1])
        ax.set_xlim(xlim[0], xlim[1])

    return fig, axs


def sketchXY( log: Logger, cases: List[Case], **kwargs ) -> Tuple[plt.Figure, List[plt.Axes]]:
    """
    """
    grid = list(_calculate_grid(cases, **kwargs))
    lcount = max([len(c.shape) for c in cases])
    grid[0] = grid[0] * lcount

    log.info(f'Generating an image grid of: {grid[0]}x{grid[1]}\n')
    fsize = (kwargs.pop('width', 7.5 * grid[1]),
             kwargs.pop('hight', 7.5 * grid[0] / lcount))
    fig = plt.figure(figsize=fsize)
    axs = []
    # ylim, xlim = [0, 0], [0, 0]
    for i, case in enumerate(cases):
        position = (int(i / grid[1]), i % grid[1])
        title = '{0}_{1:03d}'.format(case['configuration.name'], i + 1)
        log.info(f'Showing {title}-{case.architecture_str} in position: {position[0]}x{position[1]}\n')

        lcaxs = []
        for xx in range(lcount):
            p = list(copy.deepcopy(position))
            p[0] += xx
            lcaxs.append(plt.subplot2grid(grid, p, fig=fig))
        axs.extend(lcaxs)
        plot_case_sketch_vertical(case, lcaxs,
                                  kwargs.pop('connections', False),
                                  kwargs.pop('beta_fill', 'red'),
                                  kwargs.pop('beta_edge', 'black'),
                                  kwargs.pop('alpha_fill', 'blue'),
                                  kwargs.pop('alpha_edge', 'black'))
        for xx in lcaxs:
            xx.set_title(title + ' - ' + xx.get_title())

    return fig, axs


def _calculate_grid( cases: List[Case], **kwargs):
    ncases = len(cases)
    columns = kwargs.pop('columns', 2 if ncases > 1 else 1)
    return (int(ncases / columns) + ncases % columns, columns)
