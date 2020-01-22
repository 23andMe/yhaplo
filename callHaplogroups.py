#!/usr/bin/env python
#
# David Poznik
# 2016.01.08
# callHaplogroups.py
#
# yHaplo driver script
#
# To run: python -m yhaplo.callHaplogroups
#----------------------------------------------------------------------
from __future__ import absolute_import
from .config import Config
from .sample import Sample
from .tree import Tree

def callHaplogroups(useDefaultCmdLineArgs=False,
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

def callHaplogroupsOnResidList(residList):
    '''calls haplogroups over a list of 23andMe research IDs.
        returns a dictionary: key=resid, value=dictionary of results'''
    
    callHaplogroups(useDefaultCmdLineArgs=True,
                    suppressOutputAndLog=True,
                    residList=residList)

    resultsDict = dict()
    for sample in Sample.sampleList:
        resultsDict[sample.ID] = {
            'yhaplo:haplogroup': sample.haplogroup,
            'yhaplo:hgSNP':      sample.hgSNP,
            'yhaplo:hgSNPobs':   sample.hgSNPobs,
        }
        
    return resultsDict

if __name__ == '__main__':
    callHaplogroups()
