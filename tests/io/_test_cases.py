# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries

# External Libraries
import pytest
from marshmallow import ValidationError

# This Library
from topobuilder.io import CaseSchema


class TestCases( object ):
    """
    Test reading case definitions.
    """

    def test_empty( self ):
        input = {}
        error = {'configuration': ['Configuration data is required'],
                 'topology': ['A topological definition is required']}
        with pytest.raises(ValidationError) as message:
            CaseSchema().load(input)
        assert message.value.message == error

    def test_empty_config( self ):
        input = {'configuration': {}}
        error = {'configuration': {'name': ['A case identifier is required']},
                 'topology': ['A topological definition is required']}
        with pytest.raises(ValidationError) as message:
            CaseSchema().load(input)
        assert message.value.message == error
