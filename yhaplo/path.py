# David Poznik
# 2016.06.08
# path.py
#
# Defines the Path class.
#----------------------------------------------------------------------
from __future__ import absolute_import
from collections import deque

from .snp import SNP

class Path(object):
    '''
    An instance of this class represents a path through a tree.
    It stores the next node to visit, a list of SNPs observed in the derived state,
    the most derived SNP observed, and the number of ancestral alleles encountered.
    '''
    
    def __init__(self, node):
        self.node = node
        self.derSNPlist = list()
        self.mostDerivedSNP = None
        self.numAncestral = 0
        self.initPushThroughVars()
        
    def initPushThroughVars(self):
        '''
        initializes variables that track progress subsequent to pushing through
        a branch with 1 ancestral and 0 derived alleles
        '''
        
        self.nodeWhenPushedThrough = None
        self.mostDerivedSNPWhenPushedThrough = None
        self.numAncSincePushThrough = 0
        self.numDerSincePushThrough = 0
        
    def setPushThroughVars(self):
        'set memory of pushthough state to current state'

        self.nodeWhenPushedThrough = self.node
        self.mostDerivedSNPWhenPushedThrough = self.mostDerivedSNP
        
    def updatePushThroughVars(self, numAncestral, numDerived):
        'update pushthrough state with data from most recent branch assessment'
        
        self.numAncSincePushThrough += numAncestral
        self.numDerSincePushThrough += numDerived
        if self.numDerSincePushThrough > self.numAncSincePushThrough:
            self.initPushThroughVars()
        
    def copyAllAttributesOtherThanNode(self, other):
        'copies all attributes of another path, other than its node'
        
        self.derSNPlist = list(other.derSNPlist)
        self.mostDerivedSNP = other.mostDerivedSNP
        self.numAncestral = other.numAncestral
        
        self.nodeWhenPushedThrough = other.nodeWhenPushedThrough
        self.mostDerivedSNPWhenPushedThrough = other.mostDerivedSNPWhenPushedThrough
        self.numAncSincePushThrough = other.numAncSincePushThrough
        self.numDerSincePushThrough = other.numDerSincePushThrough
    
    def __str__(self):
        return '%d %d\n%s\n%s\n' % (self.numAncestral, self.numDerived,
                                    self.nodeString, self.snpString)

    # properties
    #----------------------------------------------------------------------
    @property
    def hasPushedThrough(self):
        'whether or not this path has pushed through a branch with 1 ancestral and 0 derived'
        
        return self.nodeWhenPushedThrough is not None
    
    @property
    def nodeString(self):
        'string concatenation of nodes visited'
        
        return ' '.join([node.label for node in self.node.backTracePath()])
    
    @property
    def numDerived(self):
        'number of derived SNPs in the list'

        return len(self.derSNPlist)
    
    @property
    def snpString(self):
        'string concatenation of derived SNPs observed'
        
        return ' '.join([snp.label for snp in self.derSNPlist])

    # regular methods
    #----------------------------------------------------------------------
    def betterThan(self, other):
        'evaluates whether this path is better than another'
        
        return (other is None
                or self.numDerived  > other.numDerived
                or (self.numDerived == other.numDerived
                    and self.numAncestral < other.numAncestral))

    def fork(self, nodeList):
        '''
        returns a deque of paths, each of which is identical to self
        but with a new current node
        '''
        
        pathDeque = deque()
        for node in nodeList:
            path = Path(node)
            path.copyAllAttributesOtherThanNode(self)
            pathDeque.append(path)
            
        return pathDeque

    def revertIfPushedThroughTooFar(self):
        '''
        if the path has pushed through a branch with 1 ancestral and 0 derived
        and, after doing so, it has encountered just one derived allele and a nonzero
        number of ancestral alleles, revert the path to its state before pushing through
        '''
        
        if (self.hasPushedThrough
                and self.numAncSincePushThrough  > 0
                and self.numDerSincePushThrough == 1):
            self.node = self.nodeWhenPushedThrough
            del(self.derSNPlist[-1])
            self.mostDerivedSNP = self.mostDerivedSNPWhenPushedThrough
            self.numAncestral -= self.numAncSincePushThrough
            self.initPushThroughVars()
    
    def updateWithBranchAssessment(self, ancSNPlist, derSNPlist):
        '''
        extends derived SNP list, sets most derived SNP, and adds number of
        ancestral alleles seen. also, manages tracking of whether or not path
        has pushed through an (anc,der)==(1,0) branch
        '''
            
        numAncestral, numDerived = len(ancSNPlist), len(derSNPlist)
        self.numAncestral += numAncestral
        self.derSNPlist.extend(derSNPlist)
        if derSNPlist:
            self.mostDerivedSNP = SNP.mostHighlyRankedMarkerOnList(derSNPlist)
        
        if self.hasPushedThrough:
            self.updatePushThroughVars(numAncestral, numDerived)
        elif (numAncestral, numDerived) == (1, 0):
            self.setPushThroughVars()
    
    # static methods
    #----------------------------------------------------------------------
    @staticmethod
    def createPathDeque(nodeList):
        'returns a deque of paths, with the node of each initialized from the given list'
        
        pathDeque = deque()
        for node in nodeList:
            pathDeque.append(Path(node))
            
        return pathDeque

    @staticmethod
    def postProcessPathListAndSelectBest(pathList):
        'post-processes each path in list and returns the best one'
        
        for path in pathList:
            path.revertIfPushedThroughTooFar()
            
        return Path.bestPathInList(pathList)
    
    @staticmethod
    def bestPathInList(pathList):
        'selects the best from a list of paths'
        
        bestPath = None
        for path in pathList:
            if path.betterThan(bestPath):
                bestPath = path

        return bestPath
