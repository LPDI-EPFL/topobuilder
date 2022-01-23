# -*- coding: utf-8 -*-
"""
.. codeauthor:: Zander Harteveld  <zandermilanh@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""

import os
import glob
import argparse
import numpy as np
import pandas as pd

import proq4


def cli_parser():
    """
    Create a CLI parser.
    :return: the parser object.
    """
    parser = argparse.ArgumentParser(description="Calculate lDDT using ProQ4 neural net.",
                                     epilog="Takes a PDB structure and calculates the corresponding"
                                            " lDDT score. ProQ4 must be installed.",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--prefix', '-p', type=str, nargs=1,
                        help="Prefix for output file.")
    parser.add_argument('--pdb_folder', '-f', type=str, nargs='+',
                        help="Path to folder containing the pdb files.")
    parser.add_argument('--output_folder', '-o', type=str, nargs='+',
                        help="Path to folder for generated output files.")
    #parser.add_argument('--type', '-t', type=str, nargs=1, default="naive",
    #                     help="Model type. (default: naive)")
    #parser.add_argument('--stage', '-s', type=str, nargs=1, default="design",
    #                     help="Backbone stage. (default: design)")
    return parser


def parse_args(parser):
    """
    Parse the arguments of a parser object.
    :param parser: the parser object.
    :return: the specified arguments
    """
    args = parser.parse_args()
    return args


def pdb2fasta(infile, outfile):
    """
    Creates a fasta file from a pdb file.
    """
    letters = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H',
               'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
               'TYR':'Y','VAL':'V'}

    with open(infile, "r") as infl:
        base = os.path.basename(infile).replace('.pdb', '')
        with open(outfile, "w") as otfl:
            prev = '-1'
            otfl.write(f'> {base}\n')
            for line in infl:
                toks = line.split()
                if len(toks)<1: continue
                if toks[0] != 'ATOM': continue
                if toks[2] == 'CA':
                    otfl.write('%c' % letters[toks[3]])

            otfl.write('\n')


def main():
    """
    Main entry point.
    """
    # Arguments
    args = parse_args(cli_parser())
    prefix = args.prefix[0]
    pdb_folder = args.pdb_folder[0]
    output_folder = args.output_folder[0]
    #type_ = args.type[0]
    #stage_ = args.stage[0]

    # Load model
    model = proq4.get_proq4()

    # Compute and save
    d = {#"type": [],
         #"stage": [],
         "description": [], "lDDT": [], "global lDDT": [], "classes": []}
    for pdbfile in glob.iglob(f'{pdb_folder}/*.pdb'):
        description = os.path.basename(pdbfile).replace(".pdb", "")
        #if not glob.iglob('{}/*.fasta'.format(pdb_folder)):
        fastafile = pdbfile.replace(".pdb", ".fasta")
        if not os.path.exists(fastafile):
            pdb2fasta(pdbfile, fastafile)
        m3a = proq4.process_fasta(fastafile)
        preds, classes = proq4.predict(model, m3a, pdbfile)

        d["description"].append(description)
        d["lDDT"].append( np.array(preds) )
        d["global lDDT"].append( np.mean( np.array(preds) ) )
        d["classes"].append( np.array(classes) )
        #d["type"].append(type_)
        #d["stage"].append(stage_)

    df = pd.DataFrame(d)
    df.to_csv(f'{output_folder}/_proq4.{prefix}.csv')

if __name__ == '__main__':
    main()
