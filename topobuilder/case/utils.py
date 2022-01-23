# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import string
from pathlib import Path
from typing import Optional, Tuple, Dict, List

# External Libraries
import matplotlib.pyplot as plt
from matplotlib.path import Path as pltPath
from matplotlib.transforms import Affine2D
from matplotlib.patches import PathPatch, ArrowStyle, FancyArrowPatch

# This Library
from .case import Case, layer_hights

__all__ = ['case_template', 'plot_case_sketch', 'plot_case_sketch_vertical']


def case_template( name: str,
                   architecture: Optional[str] = None,
                   topology: Optional[str] = None,
                   corrections: Optional[Dict] = None,
                   format: str = 'yaml',
                   make_absolute: bool = False
                   ) -> Tuple[Case, Path]:
    """Generate a :class:`.Case`.

    :param str name: Identifier of the case.
    :param str architecture: Definition of unconnected, unordered secondary structure.
    :param str topology: Definition of connected,ordered secondary structure. If provided, it will
        overwrite ``architecture``.
    :param dict corrections: Corrections to apply to the default case, identified by the SSE id.
    :param str format: Format of the output file (``yaml`` or ``json``).
    :param bool make_absolute: If :data:`True`, coordinates and shifts are
        defined as absolute positions.

    :return: :class:`.Case` and :class:`Path` generated filename.
    """

    # Create the case
    case = Case(name)

    # Architecture-defined case
    case = case.add_architecture(architecture)

    # Topology-defined case (architecture + connectivity)
    # If match, it is appended to the pre-defined architecture
    case = case.add_topology(topology)

    # Apply corrections if any
    case = case.apply_corrections(corrections)

    if make_absolute:
        case = case.cast_absolute()

    # Output
    outfile = case.write(name, format)

    return case, outfile


def plot_case_sketch( case: Case,
                      ax: Optional[plt.Axes] = None,
                      connections: Optional[bool] = False,
                      beta_fill: Optional[str] = 'red',
                      beta_edge: Optional[str] = 'black',
                      alpha_fill: Optional[str] = 'blue',
                      alpha_edge: Optional[str] = 'black',
                      connection_edge: Optional[str] = None
                      ) -> plt.Axes:
    """
    """
    def make_triangle(y, x, rot_deg, fcolor, ecolor, scale):
        unit_triangle = pltPath.unit_regular_polygon(3)
        path = pltPath(unit_triangle.vertices * scale, unit_triangle.codes)
        trans = Affine2D().translate(x, y).rotate_deg_around(x, y, rot_deg)
        t_path = path.transformed(trans)
        patch = PathPatch(t_path, facecolor=fcolor, edgecolor=ecolor, zorder=2)
        return patch

    if ax is None:
        fig = plt.figure()
        ax = plt.subplot2grid((1, 1), (0, 0), fig=fig)
    ax.set_aspect('equal', adjustable='box')

    shp = case.center_shape
    margin = 4
    xmax = max([shp[l]['right'] for l in shp]) + margin
    xmin = min([shp[l]['left'] for l in shp]) - margin
    ymax = 0
    ymin = 0

    for layer in case.cast_absolute()['topology.architecture']:
        for sse in layer:
            ymax = sse['coordinates']['z'] if sse['coordinates']['z'] > ymax else ymax
            ymin = sse['coordinates']['z'] if sse['coordinates']['z'] < ymin else ymin
            rotation = 180 if sse['tilt']['x'] > 90 and sse['tilt']['x'] < 270 else 0
            if sse['type'] == 'H':
                c = plt.Circle((sse['coordinates']['x'], sse['coordinates']['z']), radius=3,
                               facecolor=alpha_fill, edgecolor=alpha_edge, zorder=2)
                ax.add_artist(c)
                p = make_triangle(sse['coordinates']['z'], sse['coordinates']['x'], rotation,
                                  lighten_color(alpha_fill, 0.5), alpha_edge, 2)
                ax.add_artist(p)
            if sse['type'] == 'E':
                p = make_triangle(sse['coordinates']['z'], sse['coordinates']['x'], rotation,
                                  beta_fill, beta_edge, 2)
                ax.add_artist(p)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymax + margin, ymin - margin)
    ax.set_xlabel('X')
    ax.set_ylabel('Z')
    ax.grid(zorder=0)
    return ax


def plot_case_sketch_vertical(case: Case,
                              axs: Optional[List[plt.Axes]] = None,
                              connections: Optional[bool] = False,
                              beta_fill: Optional[str] = 'red',
                              beta_edge: Optional[str] = 'black',
                              alpha_fill: Optional[str] = 'blue',
                              alpha_edge: Optional[str] = 'black'
                              ) -> List[plt.Axes]:
    """
    """
    def make_arrow(start, end, width, tip, fcolor, ecolor):
        arrowstyle = ArrowStyle.Simple(head_length=width + tip, head_width=width + tip, tail_width=width)
        arrow = FancyArrowPatch(start, end, arrowstyle=arrowstyle)
        return PathPatch(arrow.get_path(), edgecolor=ecolor, facecolor=fcolor, zorder=2)

    layers = case.cast_absolute()['topology.architecture']
    if axs is None:
        fig = plt.figure()
        axs = []
        for x, _ in enumerate(layers):
            axs.append(plt.subplot2grid((len(layers), 1), (x, 0), fig=fig))
    if len(axs) < len(layers):
        raise ValueError('Not enough axis to plot all layers.')

    maxt, minb = 0, 0
    for i, layer in enumerate(layers):
        tops, btms = layer_hights(case, layer)
        maxt = max(tops) if max(tops) > maxt else maxt
        minb = min(btms) if min(btms) < minb else minb

    asciiU = string.ascii_uppercase
    sizes = case.center_shape
    maxw = max(sizes[l]['right'] for l in sizes)
    minw = min(sizes[l]['left'] for l in sizes)

    for i, layer in enumerate(layers):
        tops, btms = layer_hights(case, layer)
        for isse, sse in enumerate(layer):
            rotation = True if sse['tilt']['x'] > 90 and sse['tilt']['x'] < 270 else False
            scale = 2 if sse['type'] == 'E' else 3
            tip = (scale / 2) if sse['type'] == 'E' else 0
            ecolor = beta_edge if sse['type'] == 'E' else alpha_edge
            fcolor = beta_fill if sse['type'] == 'E' else alpha_fill
            if not rotation:
                axs[i].add_artist(make_arrow([sse['coordinates']['x'], btms[isse]],
                                             [sse['coordinates']['x'], tops[isse]],
                                             scale, tip, fcolor, ecolor))
            else:
                axs[i].add_artist(make_arrow([sse['coordinates']['x'], tops[isse]],
                                             [sse['coordinates']['x'], btms[isse]],
                                             scale, tip, fcolor, ecolor))
            axs[i].set_ylim(minb, maxt)
            axs[i].set_xlim(minw - 1.7, maxw + 1.7)
            axs[i].set_title('Layer {}'.format(asciiU[i]))
            axs[i].grid(zorder=0)
            axs[i].set_xlabel('X')
            axs[i].set_ylabel('Y')
    return axs


def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    https://stackoverflow.com/a/49601444/2806632

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except Exception:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])
