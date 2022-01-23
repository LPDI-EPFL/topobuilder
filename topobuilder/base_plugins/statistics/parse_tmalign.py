# -*- coding: utf-8 -*-
"""
.. codeauthor:: Zander Harteveld  <zandermilanh@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""

import os
import sys
import argparse
import glob

import pandas as pd


def cli_parser():
    """
    Create a CLI parser.
    :return: the parser object.
    """
    parser = argparse.ArgumentParser(description="Parse TMalign files.",
                                     epilog="Reads a TMalgin file and stores the necessary"
                                            " information into a pandas dataframe.",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--tmalign_folder', '-f', type=str, nargs=1, help="the TMalign folder.")
    parser.add_argument('--out_file', '-o', type=str, nargs=1, help="the processed output file.")
    return parser


def parse_args(parser):
    """
    Parse the arguments of a parser object.
    :param parser: the parser object.
    :return: the specified arguments
    """
    args = parser.parse_args()
    return args


def read_tmalign(folder):
    """
    Read TMalign file.
    :param filename: the path to the TMalign output file.
    :return: pandas DataFrame with info.
    """
    # Set frame
    d = {'description': [], 'target': [],
         'l_description': [], 'l_target': [],
         'tm_score': [], 'rmsd': [], 'seqid': [], 'n_aligned': [],
         'alignment': []}

    # Read files
    for filename in glob.iglob(f'{folder}_trRosetta*'):
       with open(filename, 'r') as f:
          lines = f.readlines()

       # Find
       alignment_str = ''
       tmscore, _alignment = None, 0.
       for line in lines:
          # Note that chain 1 gets superimposed onto chain 2,
          # thus we save chain 1 under name2 and viceversa
          if line.startswith('Name of Chain_1'):
             name2 = [s.strip() for s in line.split() if s.endswith('.pdb')][0]
             #name2 = [s.strip() for s in line.split()][-1]
             name2 = os.path.basename(name2).replace('.pdb', '')
          if line.startswith('Name of Chain_2'):
             name1 = [s.strip() for s in line.split() if s.endswith('.pdb')][0]
             #name1 = [s.strip() for s in line.split()][-1]
             name1 = os.path.basename(name1).replace('.pdb', '')
          if line.startswith('Length of Chain_1'):
             l2 = int(line.split()[-2].strip())
          if line.startswith('Length of Chain_2'):
             l1 = int(line.split()[-2].strip())
          if line.startswith('Aligned'):
             infos = line.split('=')
             al    = int(infos[1].split(',')[0].strip())
             rmsd  = float(infos[2].split(',')[0].strip())
             seqid = float(infos[-1].strip())
          if not tmscore and line.startswith('TM-score'):
             tmscore = float(line.split()[1].strip())
          if tmscore and line.startswith('TM-score'):
             tmscore2 = float(line.split()[1].strip())
             if tmscore < tmscore2:
                tmscore = tmscore2
          if line.startswith('(":"'):
             _alignment += 1
             continue
          if _alignment > 0:
             if _alignment < 4:
                alignment_str += line.strip()
                alignment_str += '\n'
                _alignment += 1

       # Save
       d['description'].append(name2)
       d['target'].append(name1)
       d['l_description'].append(l1)
       d['l_target'].append(l2)
       d['tm_score'].append(tmscore)
       d['rmsd'].append(rmsd)
       d['seqid'].append(seqid)
       d['n_aligned'].append(al)
       d['alignment'].append(alignment_str)

    return pd.DataFrame(d)


def main():
    """
    Main entry point.
    """
    # Arguments
    args = parse_args(cli_parser())
    tmalign_folder = args.tmalign_folder[0]
    out_file = args.out_file[0]

    # Fetch and write
    #sys.stdout.write(f'Reading from: {tmalign_folder}\n')
    df = read_tmalign(tmalign_folder)
    #sys.stdout.write(f'Writing file: {out_file}.csv\n')
    df.to_csv(f'{out_file}')


if __name__ == '__main__':
    main()
