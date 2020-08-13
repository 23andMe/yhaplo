# David Poznik
# 2015.12.29
# node.py
#
# Defines the Node class.
#----------------------------------------------------------------------
from __future__ import absolute_import
from collections import deque
from operator import attrgetter
from six.moves import range

from . import utils
from .page import Page
from .snp import SNP

class Node(object):
    '''
    A node knows its:
        - parent (self.parent is None == self.isRoot())
        - depth
        - children
        - diagnostic SNPs
    
    Throughout this code, each node represents the branch that leads to it.
    '''

    tree = None
    config = None
    args = None
    errAndLog = None
    pageList = list()
    pageDict = dict()
    hgSNPset = set()

    def __init__(self, parent, tree=None):
        self.parent = parent
        if self.isRoot():
            Node.setTreeConfigAndArgs(tree)
            self.depth = 0
        else:
            parent.addChild(self)
            self.depth = parent.depth + 1
            if self.depth > Node.tree.maxDepth:
                Node.tree.maxDepth = self.depth

        self.haplogroup = ''            # see setLabel  | ycc haplogroup name     e.g., R1b1c
        self.label = ''                 # see setLabel  | ycc including alt names e.g., P/K2b2
        self.hgTrunc = ''               # see setLabel  | truncated haplogroup    e.g., R1b1c -> R
        self.hgSNP = ''                 # see prioritySortSNPlistAndSetHgSNP |    e.g., R-V88
        self.childList = list()         # see addChild, bifurcate, serialSplit
        self.snpList = list()           # see addSNP
        self.droppedMarkerList = list() # see addDroppedMarker
        self.page = None                # see addSNP
        self.branchLength = None        # see setBranchLength
        self.DFSrank = 0                # see setDFSrank

        if Node.args.writeContentMappings and self.isRoot():
            self.page = Node.pageDict[Node.config.rootHaplogroup]
            self.page.setNode(self)
        
    def __str__(self):
        return self.strSimple()
    
    def strSimple(self):
        'string representation: label and representative SNP'
        
        return '%-25s %s' % (self.label, self.hgSNP)
    
    def strSNPlist(self):
        'string representation: label and list of snps'
        
        snpString = ' '.join(snp.label for snp in self.snpList)
        return '%-25s %s' % (self.label, snpString)
        
    def strDotPipeDepth(self):
        'string representation: indicates depth with a series of dots and pipes'
        
        dotList = list('.' * (self.depth))
        for i in range(0, len(dotList), 5):
            dotList[i] = '|'
        return '%s%s %s' % (''.join(dotList), self.label, self.hgSNP)

    def strTreeTableRow(self):
        'string representation: one row of tree table'
        
        yccLabel = self.haplogroup
        
        if self.isRoot():
            parentDFSrank = parentHgSNP = 'root'
        else:
            parentDFSrank = str(self.parent.DFSrank)
            parentHgSNP = self.parent.hgSNP
            
        return '\t'.join([str(self.DFSrank), yccLabel, self.hgSNP, parentDFSrank, parentHgSNP])

    @property
    def mostHighlyRankedSNP(self):
        'the most highly ranked SNP'
        
        return SNP.mostHighlyRankedMarkerOnList(self.snpList)
    
    @property
    def mostHighlyRankedDroppedMarker(self):
        'the most highly ranked dropped marker'
        
        return SNP.mostHighlyRankedMarkerOnList(self.droppedMarkerList)

        
    # static methods, including class variable setters
    #----------------------------------------------------------------------
    @staticmethod
    def setTreeConfigAndArgs(tree):
        'enables Node class to know about the tree instance, config, and args'
        
        Node.tree = tree
        Node.config = tree.config
        Node.args = tree.args
        Node.errAndLog = tree.config.errAndLog
        if Node.args.writeContentMappings:
            Node.buildPageDict()

    @staticmethod
    def buildPageDict():
        '''
        builds a dictionary of 23andMe content pages. pagesFN comes from these two gdocs:
        - https://docs.google.com/spreadsheets/d/1mf86slweZEKUd5hzG2GmKGTGIpHuDipJz2u221y2zVE/edit?ts=568eb997#gid=0
        - https://docs.google.com/spreadsheets/d/1oo0sRmYFNeWikuOxcb_1obOoO35wQccmOzyGRmqDMtc/edit?ts=578578d0#gid=362797346
        '''
        
        utils.checkFileExistence(Node.config.pagesFN, 'Content pages')
        with open(Node.config.pagesFN, 'r') as pagesFile:
            pagesFile.readline()    # header
            for line in pagesFile:
                yccOld, snpName = line.strip().split()
                page = Page(yccOld, snpName)
                Node.pageList.append(page)
                
                if yccOld == Node.config.rootHaplogroup:
                    Node.pageDict[Node.config.rootHaplogroup] = page
                elif snpName != '.':
                    Node.pageDict[snpName] = page

    @staticmethod
    def truncateHaplogroupLabel(haplogroup):
        'returns first 2-5 characters of specified haplogroups and first letter of others'
        
        for numChars in range(Node.config.multiCharHgTruncMaxLen, 1, -1):
            if haplogroup[:numChars] in Node.config.multiCharHgTruncSet:
                return haplogroup[:numChars]
        
        return haplogroup[0]

    # setters, mutaters
    #----------------------------------------------------------------------
    def setLabel(self, label):
        'sets label, haplogroup, and hgTrunc'
        
        self.label = label
        labelList = label.split('/')
        
        if self.isRoot():
            self.haplogroup = self.hgTrunc = self.config.rootHaplogroup
            Node.tree.hg2nodeDict[self.haplogroup] = self
        else:
            self.haplogroup = labelList[0]
            self.hgTrunc = Node.truncateHaplogroupLabel(self.haplogroup)
        
        for key in labelList:
            Node.tree.hg2nodeDict[key] = self
            
    def setBranchLength(self, branchLength):
        'sets the branch length'
        
        self.branchLength = branchLength
        
    def setDFSrank(self, DFSrank):
        'set depth-first search rank'
        
        self.DFSrank = DFSrank
        
    def addSNP(self, snp):
        'appends a snp to the snp list'
        
        self.snpList.append(snp)
        if snp.label in Node.pageDict:
            self.page = Node.pageDict[snp.label]
            self.page.setNode(self)
    
    def addDroppedMarker(self, droppedMarker):
        'appends a dropped marker to the list'
        
        self.droppedMarkerList.append(droppedMarker)
    
    def prioritySortSNPlistAndSetHgSNP(self):
        '''
        first, sorts snp list (or dropped marker list) by priority ranking.
        then, sets reresentative-SNP-based label: self.hgSNP
        the standard form incudes the truncated haplogroup label
        and the label of a representative SNP, separated by a hyphen (e.g. R-V88).
        '''
        
        # root: no markers
        if self.isRoot():
            self.hgSNP = self.haplogroup
            
        # normal case
        elif self.snpList:
            self.snpList = SNP.prioritySortMarkerList(self.snpList)
            self.hgSNP = self.mostHighlyRankedSNP.hgSNP
            
        # backup: use discared marker name
        elif self.droppedMarkerList:
            self.droppedMarkerList = SNP.prioritySortMarkerList(self.droppedMarkerList)
            markerName = self.mostHighlyRankedDroppedMarker.name
            self.hgSNP = '%s-%s' % (self.hgTrunc, markerName)
            
        # no markers to use
        else:
            if self.parent.hgSNP:
                symbol = '*' if self.isLeaf() else '+'
                self.hgSNP = self.parent.hgSNP + symbol
                
                # uniquify if necessary
                if self.hgSNP in Node.hgSNPset:
                    i = 1
                    hgSNPuniqe = '%s%d' % (self.hgSNP, i)
                    while hgSNPuniqe in Node.hgSNPset:
                        i += 1
                        hgSNPuniqe = '%s%d' % (self.hgSNP, i)
                    
                    self.hgSNP = hgSNPuniqe
            else:
                Node.errAndLog('WARNING. Attempted to set star label, ' +
                               'but parent.hgSNP not set yet: %s\n' % self.haplogroup)
                self.hgSNP = self.haplogroup
                
        Node.hgSNPset.add(self.hgSNP)
        
    # queries
    #----------------------------------------------------------------------
    def isRoot(self):
        return self.parent is None

    def isLeaf(self):
        return len(self.childList) == 0
    
    def getBranchLength(self, alignTips=False, platformVersion=None):
        if self.branchLength:
            return self.branchLength
        elif alignTips and self.isLeaf():
            return Node.tree.maxDepth - self.depth + 1
        elif alignTips:
            return 1
        elif platformVersion:
            branchLength = 0
            for snp in self.snpList:
                if snp.isOnPlatform(platformVersion):
                    branchLength += 1
            return branchLength
        else:
            return None
    
    def backTracePath(self):
        'returns a list of nodes from root to self'
        
        nodeList = [self]
        parent = self.parent
        while parent is not None:
            nodeList.append(parent)
            parent = parent.parent
        nodeList.reverse()
        return nodeList

    def assessGenotypes(self, sample):
        '''
        assess an individual's genotypes with respect to self.snpList
        returns two lists of snps. those for which:
            - ancestral genotypes were observed
            - derived genotypes were observed
        '''
        
        genotypedSnpList = [snp for snp in self.snpList
                            if snp.position in sample.pos2genoDict]
        ancSNPlist, derSNPlist = list(), list()
        listAllGenotypes = Node.args.haplogroupToListGenotypesFor == self.haplogroup

        for snp in genotypedSnpList:
            geno = sample.pos2genoDict[snp.position]
            
            if snp.isAncestral(geno):
                ancSNPlist.append(snp)
            elif snp.isDerived(geno):
                derSNPlist.append(snp)

            if listAllGenotypes:
                derivedFlag = '*' if snp.isDerived(geno) else ''
                Node.config.hgGenosFile.write('%-8s %s %s %s\n' %
                                              (sample.ID, snp, geno, derivedFlag))
        
        return ancSNPlist, derSNPlist

    # children
    #----------------------------------------------------------------------
    def addChild(self, child):
        'appends a child to the child list'
        
        self.childList.append(child)
    
    def serialSplit(self, targetHaplogroup):
        'serially split node until there is a spot for the target haplogroup'
        
        currentNode = self
        startLength = len(self.haplogroup)
        endLength = len(targetHaplogroup)
        for strLen in range(startLength, endLength):
            nextNode = None
            targetHgSubstring = targetHaplogroup[:(strLen+1)]
            if currentNode.numChildren < 2:
                currentNode.bifurcate()
            for node in currentNode.childList:
                if node.haplogroup == targetHgSubstring:
                    nextNode = node
            if nextNode is None:
                nextNode = Node(parent = currentNode)
                nextNode.setLabel(targetHgSubstring)
                currentNode.sortChildren()

            currentNode = nextNode

        return currentNode

    @property
    def numChildren(self):
        return len(self.childList)
    
    def bifurcate(self):
        'split a node and return the two children'
        
        leftChild = Node(parent = self)
        rightChild = Node(parent = self)
        if self.haplogroup[-1].isalpha():
            leftChild.setLabel(self.haplogroup + '1')
            rightChild.setLabel(self.haplogroup + '2')
        else:
            leftChild.setLabel(self.haplogroup + 'a')
            rightChild.setLabel(self.haplogroup + 'b')
        return leftChild, rightChild
    
    def sortChildren(self):
        self.childList = sorted(self.childList, key=attrgetter('haplogroup'))
    
    def reverseChildren(self):
        self.childList.reverse()
        
    # tree traversals
    #----------------------------------------------------------------------
    def writeBreadthFirstTraversal(self, bfTreeFile):
        'writes breadth-first traversal'
        
        bfTreeFile.write('%s\n' % self.strDotPipeDepth())
        nodeDeque = deque(self.childList)
        while nodeDeque:
            node = nodeDeque.popleft()
            bfTreeFile.write('%s\n' % node.strDotPipeDepth())
            nodeDeque.extend(node.childList)
        
    def getDepthFirstNodeList(self):
        'wrapper function for recursive depth-first pre-order traversal'
        
        depthFirstNodeList = [self]
        self.traverseDepthFirstPreOrderRecursive(depthFirstNodeList)
        return depthFirstNodeList
        
    def traverseDepthFirstPreOrderRecursive(self, depthFirstNodeList):
        'recursively appends each node in depth-first pre order'
        
        for child in self.childList:
            depthFirstNodeList.append(child)
            child.traverseDepthFirstPreOrderRecursive(depthFirstNodeList)
    
    def mrca(self, otherNode):
        'returns the most recent common ancestor of this node and another'
        
        if self.depth < otherNode.depth:
            higherNode, lowerNode = self, otherNode
        else:
            higherNode, lowerNode = otherNode, self
        
        while higherNode.depth < lowerNode.depth:
            lowerNode = lowerNode.parent
        while lowerNode != higherNode:
            lowerNode = lowerNode.parent
            higherNode = higherNode.parent
        
        return higherNode

    # writing tree to file in Newick format
    #----------------------------------------------------------------------
    def writeNewick(self, newickFN,
                    useHgSNPlabel=False, alignTips=False, platformVersion=None):
        'write Newick string for the subtree rooted at this node'

        if not Node.config.suppressOutputAndLog:
            with open(newickFN, 'w') as outFile:
                outFile.write('%s;\n' %
                    self.buildNewickStringRecursive(
                        useHgSNPlabel, alignTips, platformVersion))
            
            if alignTips:
                treeDescriptor = 'aligned '
            elif platformVersion:
                treeDescriptor = 'platform v%d ' % platformVersion
            else:
                treeDescriptor = ''
                
            if useHgSNPlabel:
                labelType = 'representative-SNP'
            else:
                labelType = 'YCC'
                
            Node.errAndLog('Wrote %stree with %s labels:\n    %s\n\n' %
                           (treeDescriptor, labelType, newickFN))

    def buildNewickStringRecursive(self,
            useHgSNPlabel=False, alignTips=False, platformVersion=None):
        'recursively builds Newick string for the subtree rooted at this node'
        
        if not self.isLeaf():
            childStringList = list()
            for child in self.childList[::-1]:
                childString = child.buildNewickStringRecursive(useHgSNPlabel,
                                                               alignTips, platformVersion)
                childStringList.append(childString)
            treeStringPart1 = '(%s)' % ','.join(childStringList)
        else:
            treeStringPart1 = ''

        branchLabel = self.hgSNP if useHgSNPlabel else self.label
        branchLength = self.getBranchLength(alignTips, platformVersion)
        if alignTips:
            branchString = '%s:%d' % (branchLabel, branchLength)
        elif branchLength is None or (self.isLeaf() and branchLength == 0):
            branchString = branchLabel
        elif branchLength > 0:
            branchString = '%s|%d:%d' % (branchLabel, branchLength, branchLength)
        else:
            branchString = ':0.5'
        
        treeString = '%s%s' % (treeStringPart1, branchString)
        return treeString
