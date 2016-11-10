#!/usr/bin/env python
#
# David Poznik
# 2016.01.08
# callHaplogroups.py
#
# Driver script: builds tree and calls haplogroups for individuals
#----------------------------------------------------------------------
from config import Config
from tree import Tree
from sample import Sample

config = Config()
tree = Tree(config)
Sample.callHaplogroups(config, tree)
