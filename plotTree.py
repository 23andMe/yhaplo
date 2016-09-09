#!/usr/bin/env python
#
# David Poznik
# 2016.1.7
# plotTree.py
#----------------------------------------------------------------------
import argparse
import sys

try:
    from Bio import Phylo
except ImportError: 
    sys.exit('\nERROR. Please install Biopython with the following command:\n\n\t' +
             'pip install biopython\n')

from config import Config


DESCRIPTION = 'plots a newick tree'

parser = argparse.ArgumentParser(description=DESCRIPTION, 
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-d', '--draw', action='store_true', default=False,
                    help = 'draw tree, rather than printing ascii version')
parser.add_argument('-n', '--newickFN', 
                    type=str, default=Config.primaryTreeFN,
                    help = 'name of file containing newick tree to plot')

args = parser.parse_args()

#----------------------------------------------------------------------
phyloTree = Phylo.read(args.newickFN, 'newick')
if args.draw:
    Phylo.draw(phyloTree)
else:
    Phylo.draw_ascii(phyloTree)
