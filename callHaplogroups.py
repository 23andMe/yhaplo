#!/usr/bin/env python
#
# David Poznik
# 2016.01.08
# callHaplogroups.py
#
# Main engine: builds tree and calls haplogroups for individuals.
#----------------------------------------------------------------------
from config import Config
from tree import Tree
from sample import Sample

description = '''
This software identifies the Y-chromosome haplogroup of each male in a sample of 
one to millions. Sequence data will yield the most highly resolved classifications, 
but the algorithm also works well with chip-based genotype data, provided a reasonable 
number of phylogenetically informative sites have been assayed. 
'''

config = Config(description)
tree   = Tree(config)
Sample.callHaplogroups(config, tree)
