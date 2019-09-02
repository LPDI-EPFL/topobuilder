TopoBuilder
===========

*TopoBuilder* is a **Python2 only** library for the construction of tailored protein scaffolds around
functional structural motifs of interest.

Provided a PDB-formatted file from the **Protein Data Bank** (PDB_) with the functional motif of interest
and a **JSON** configuration file with the expected topology, *TopoBuilder* will generate all the necessary
files to execute the FunFolDes_ protocol of the Rosetta_ suite in a non-mpi, SLURM-managed cluster.

Installation
------------

*TopoBuilder* can be installed for **Python 2.7** from the git repository through pip:

.. code-block:: bash

  pip install https://github.com/lpdi-epfl/topobuilder/archive/v1.0.zip


Execute
-------

To see all the options to execute *TopoBuilder*, write:

.. code-block:: bash

  python -m topobuilder -h


*TopoBuilder* builds the topologies by initially placing secondary structures in a X-Z grid,
where Z represents depth and X represents width. The Y axis (which is high) can also be specified,
but it is not necessary. Initially, the protocol, given a position and number of secondary structures
(architecture), will try all posible connectivities (topologies) fitting that architecture as long as
the motif allows it (multi-segment motifs, by their nature, limit combinatorial possibilities).
The best way to explain how to work with *TopoBuilder* is through an example.
Let's say on wants to build a topology carrying the 4e10 epitope from HIV (2FX7_).

To simplify, we will try a **4-helix bundle**. For that, we would write a **JSON** input such as:

.. code-block:: json

  {
    "config": {
      "name": "4hb"
    },
    "layers" : [
      [
        { "type" : "H",
          "length" : 16,
        },
        { "type" : "H",
          "length" : 16,
        }
      ],
      [
        { "type" : "H",
          "ref": "2fx7.helix"
        },
        { "type" : "H",
          "length" : 16,
        }
      ]
    ],
    "motifs": [
      {
        "id": "2fx7",
        "pdbfile": "2fx7.pdb",
        "chain": "P",
        "segments": [
          {
            "ini": 671,
            "end": 686,
            "id": "helix"
          }
        ]
      }
    ]
  }


``config``, ``layers`` and ``motifs`` are the top, mandatory fields.

``config``
**********

The only mandatory parameter here is ``name``, which identifies the full execution.
Other parameters that can be provided but have default values are:

* ``default_z``: Default depth between secondary structure layers. (default=11)
* ``default_x_h``: Default width between helices in the same layer. (default=11)
* ``default_x_e``: Default width between beta strands in the same layer. (default=5)
* ``link_dist``: Defalut distance between secondary structure to consider connecting them.
* ``connectivity``: If provided, create a given connectivity instead of trying all possible.
  Connectivity should be defined as a string in FORM_ format, in which each secondary structure
  is defined by `<layer_id><layer_position><SSE_type>`; where `<layer_id>` is an uppercase letter
  starting in A, `<layer_position>` is an integer starting in 1 and `<SSE_type>` is either (H) helix
  or (E) beta. We will see how this looks like in the results from the example execution.
* ``l_linkers``: If provided as a list of numbers with length=number of structures - 1, it will
  setup those as the loop lengths, otherwise the protocol will calculate the most likely lengths for the loops.


.. _PDB: https://www.rcsb.org/
.. _FunFolDes: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006623
.. _Rosetta: https://www.rosettacommons.org/
.. _2FX7: https://www.rcsb.org/structure/2FX7
.. _FORM: https://www.sciencedirect.com/science/article/pii/S0969212609002950