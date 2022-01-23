.. _utils:

.. currentmodule:: topobuilder

Utility Functions
=================

Structure
---------

.. autofunction:: topobuilder.utils.build_pdb_object
.. autofunction:: topobuilder.utils.get_loop_length
.. autofunction:: topobuilder.utils.pdb_geometry_from_rules
.. autofunction:: topobuilder.utils.make_angles_and_distances
.. autofunction:: topobuilder.utils.pick_motif


Plot
----

.. autofunction:: topobuilder.utils.plot_fragment_templates
.. autofunction:: topobuilder.utils.plot_loop_length_distribution
.. autofunction:: topobuilder.utils.plot_match_bin
.. autofunction:: topobuilder.utils.plot_geometric_distributions
.. autofunction:: topobuilder.utils.plot_angle_network


MASTER
------

.. autofunction:: topobuilder.utils.pds_database
.. autofunction:: topobuilder.utils.createPDS
.. autofunction:: topobuilder.utils.master_best_each
.. autofunction:: topobuilder.utils.master_fixedgap
.. autofunction:: topobuilder.utils.master_groupedgap
.. autofunction:: topobuilder.utils.master_nogap


Plugins
-------

.. autofunction:: topobuilder.utils.plugin_title
.. autofunction:: topobuilder.utils.plugin_case_count
.. autofunction:: topobuilder.utils.plugin_warning
.. autofunction:: topobuilder.utils.plugin_filemaker
.. autofunction:: topobuilder.utils.plugin_imagemaker
.. autofunction:: topobuilder.utils.plugin_filereader
.. autofunction:: topobuilder.utils.plugin_bash
.. autofunction:: topobuilder.utils.plugin_conditions


Rosetta
-------

.. autofunction:: topobuilder.utils.constraint_minimization
.. autofunction:: topobuilder.utils.constraint_design
.. autofunction:: topobuilder.utils.funfoldes
.. autofunction:: topobuilder.utils.hybridize
.. autofunction:: topobuilder.utils.rosettascript
.. autofunction:: topobuilder.utils.SELECTOR_SecondaryStructure
.. autofunction:: topobuilder.utils.MOVER_SetSecStructEnergies
.. autofunction:: topobuilder.utils.MOVER_PeptideStubMover
.. autofunction:: topobuilder.utils.PROTOCOL_LayerDesign
.. autofunction:: topobuilder.utils.PROTOCOL_BasicFilters


Slurm
-----

.. autofunction:: topobuilder.utils.submit_slurm
.. autofunction:: topobuilder.utils.submit_nowait_slurm
.. autofunction:: topobuilder.utils.control_slurm_file
