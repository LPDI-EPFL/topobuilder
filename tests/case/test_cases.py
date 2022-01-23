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

# External Libraries
import pytest
from marshmallow import ValidationError
import matplotlib as mpl
if os.environ.get('DISPLAY', '') == '':
    mpl.use('Agg')
import matplotlib.pyplot as plt

# This Library
from topobuilder.case import Case, plot_case_sketch


class TestCases( object ):
    """
    Test reading case definitions.
    """
    def setup( self ):
        self.datadir = Path.cwd().joinpath('..', 'data').resolve()

    def test_empty( self ):
        input = {}
        error = {'configuration': ['Configuration data is required'],
                 'topology': ['A topological definition is required']}
        with pytest.raises(ValidationError) as message:
            Case(input)
        assert message.value.messages == error

    def test_empty_config( self ):
        input = {'configuration': {}}
        error = {'configuration': {'name': ['A case identifier is required']},
                 'topology': ['A topological definition is required']}
        with pytest.raises(ValidationError) as message:
            Case(input)
        assert message.value.messages == error

    def test_minimal( self ):
        c = Case('test_minimal')
        assert c['configuration.name'] == 'test_minimal'
        assert c['topology.architecture'] == [[]]
        assert 'architecture' not in c
        assert c.shape == ()
        c.write(self.datadir)
        c.write(self.datadir, format='json')

    @pytest.mark.mpl_image_compare(baseline_dir=os.path.join('..', 'baseline_images'),
                                   filename='architecture.png')
    def test_architecture( self ):
        c = Case('test_architecture')
        c = c.add_architecture('5E.2H')
        assert c.shape == (5, 2)
        assert c.shape_len == ((7, 7, 7, 7, 7), (13, 13))
        assert c.architecture_str == '5E.2H'
        assert c.center_shape == {'A': {'bottom': 0, 'hight': 0, 'left': 0,
                                        'right': 19.4, 'top': 0, 'width': 19.4},
                                  'B': {'bottom': 0, 'hight': 0, 'left': 0,
                                        'right': 10.0, 'top': 0, 'width': 10.0}}
        with pytest.raises(ValidationError) as message:
            c = c.add_architecture('5E:8:8:7:7:7.2H:18:19')
        assert message.value.messages == ['An arquitecture is already defined.']

        c = Case('test_architecture')
        c = c.add_architecture('5E:8:8:7:7:7.2H:18:19')
        assert c.shape == (5, 2)
        assert c.shape_len == ((8, 8, 7, 7, 7), (18, 19))
        assert c.connectivity_count == 0
        assert c.connectivities_str == []
        fig = plt.figure(figsize=(15, 5))
        ax1 = plt.subplot2grid((1, 3), (0, 0), fig=fig)
        plot_case_sketch(c, ax1)
        ax2 = plt.subplot2grid((1, 3), (0, 1), fig=fig)
        plot_case_sketch(c.apply_corrections({'B': {'xalign': 'center'}}), ax2)
        ax3 = plt.subplot2grid((1, 3), (0, 2), fig=fig)
        plot_case_sketch(c.apply_corrections({'B': {'xalign': 'right'}}), ax3)
        return fig

    @pytest.mark.mpl_image_compare(baseline_dir=os.path.join('..', 'baseline_images'),
                                   filename='topology.png')
    def test_topology( self ):
        c = Case('test_topology')
        c = c.add_topology('A2E.A1E.B1H.A3E.B2H.A5E.A4E')
        assert c.shape == (5, 2)
        assert c.shape_len == ((7, 7, 7, 7, 7), (13, 13))
        assert c.architecture_str == '5E.2H'
        assert c.connectivity_count == 1
        assert c.connectivities_str == ('A2E.A1E.B1H.A3E.B2H.A5E.A4E', )
        c = c.add_topology('A2E.A1E.B1H.A3E.B2H.A5E.A4E')
        with pytest.raises(ValidationError) as message:
            c = c.add_topology('A2E.A1E.B1H.A3E.B2H.A5E.A4E.B3H')
        assert message.value.messages == ['Provided topology does not match existing architecture.']
        assert c.connectivity_count == 1
        assert c.connectivities_str == ('A2E.A1E.B1H.A3E.B2H.A5E.A4E', )
        cs = c.apply_topologies()
        assert len(cs) == 1
        assert cs[0].shape_len == c.shape_len
        assert cs[0].shape == c.shape
        assert cs[0].connectivities_str == c.connectivities_str

        c = Case('test_topology')
        c = c.add_architecture('5E:8:8:7:7:7.2H:18:19')
        c = c.add_topology('A2E8.A1E8.B1H18.A3E7.B2H19.A5E7.A4E7')
        c = c.add_topology('A2E8.A1E8.B1H18.A3E7.B2H19.A4E7.A5E7')
        c = c.add_topology('A1E8.A2E8.B1H18.A3E7.B2H19.A4E7.A5E7')
        assert c.shape == (5, 2)
        assert c.shape_len == ((8, 8, 7, 7, 7), (18, 19))
        assert c.architecture_str == '5E.2H'
        assert c.connectivity_count == 3
        assert c.connectivities_str == ('A2E.A1E.B1H.A3E.B2H.A5E.A4E',
                                        'A2E.A1E.B1H.A3E.B2H.A4E.A5E',
                                        'A1E.A2E.B1H.A3E.B2H.A4E.A5E')
        cs = c.apply_topologies()
        assert len(cs) == 3
        assert cs[0].connectivity_count == 1
        assert cs[0].connectivities_str == ('A2E.A1E.B1H.A3E.B2H.A5E.A4E', )
        assert cs[1].connectivity_count == 1
        assert cs[1].connectivities_str == ('A2E.A1E.B1H.A3E.B2H.A4E.A5E', )
        assert cs[2].connectivity_count == 1
        assert cs[2].connectivities_str == ('A1E.A2E.B1H.A3E.B2H.A4E.A5E', )

        fig = plt.figure(figsize=(15, 5))
        ax1 = plt.subplot2grid((1, 3), (0, 0), fig=fig)
        plot_case_sketch(cs[0], ax1)
        ax2 = plt.subplot2grid((1, 3), (0, 1), fig=fig)
        plot_case_sketch(cs[1], ax2)
        ax3 = plt.subplot2grid((1, 3), (0, 2), fig=fig)
        plot_case_sketch(cs[2], ax3)
        return fig
