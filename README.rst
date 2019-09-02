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
```
pip install https://github.com/lpdi-epfl/topobuilder/archive/v1.0.zip
```

Execute
-------

To see all the options to execute *TopoBuilder*, write:
```
python -m topobuilder -h
```




.. _PDB: https://www.rcsb.org/
.. _FunFolDes: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006623
.. _Rosetta: https://www.rosettacommons.org/