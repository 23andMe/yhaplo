#!/usr/bin/env python
#
# David Poznik
# 2016.01.08
# call_haplogroups.py
#
# yhaplo driver script
#
# To run: python -m yhaplo.call_haplogroups
#----------------------------------------------------------------------
from __future__ import absolute_import

import six

from .config import Config
from .sample import Sample
from .tree import Tree

def call_haplogroups(useDefaultCmdLineArgs=False,
                     suppressOutputAndLog=False,
                     outDir=None,
                     residList=None):
    'configures run, builds tree, and calls haplogroups'
    
    config = Config(useDefaultCmdLineArgs=useDefaultCmdLineArgs,
                    suppressOutputAndLog=suppressOutputAndLog,
                    outDir=outDir,
                    residList=residList)
    tree = Tree(config)
    Sample.callHaplogroups(config, tree)

def call_haplogroups_for_resid_list(resid_list):
    '''
    calls haplogroups over a list of 23andMe research IDs.
    returns a dictionary: key=resid, value=dictionary of results
    '''
    
    call_haplogroups(useDefaultCmdLineArgs=True,
                     suppressOutputAndLog=True,
                     residList=resid_list)
    
    results_dict = dict()
    for sample in Sample.sampleList:
        results_dict[sample.ID] = {
            'yhaplo:haplogroup': six.ensure_str(sample.haplogroup),
            'yhaplo:hgSNP': six.ensure_str(sample.hgSNP),
            'yhaplo:hgSNPobs': six.ensure_str(sample.hgSNPobs),
        }
    
    return results_dict

if __name__ == '__main__':
    call_haplogroups()
