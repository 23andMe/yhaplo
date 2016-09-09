# David Poznik
# 2015.12.29
# node.py
#
# Defines the Node class.
#----------------------------------------------------------------------
import sys
from collections import deque
from operator import attrgetter
from page import Page

class Node(object):
    '''
    A node knows its:
        - parent ( self.parent is None == self.isRoot() )
        - depth
        - children
        - diagnostic SNPs
    
    Throughout this code, each node represents the branch that leads to it. 
    '''

    tree      = None
    config    = None
    args      = None
    errAndLog = None
    pageDict  = dict()

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

        self.label = ''           # see setLabel
        self.haplogroup = ''      # see setLabel
        self.hgShort = ''         # see setHgShort
        self.branchLength = None  # see setBranchLength
        self.DFSrank = 0          # see setDFSrank
        self.snpList = list()     # see addSNP
        self.childList = list()   # see addChild, bifurcate, serialSplit
        self.page = None          # see setLabel, addSNP
        
    def __str__(self):
        'string representation: indicates depth for traversal output'
        
        dotList = list('.' * (self.depth))
        for i in xrange(0, len(dotList), 5):
            dotList[i] = '|'
        return '%s%s %s' % (''.join(dotList), self.label, self.strHgShortWithSNP())
    
    def strHgShortWithSNP(self):
        'string representation: haplogroup short form with SNP. e.g., R1b-V88'
        
        if self.snpList:
            return '%s-%s' % (self.hgShort, self.mostHighlyRankedSNP.label)
        else:
            return ''
    
    def strSNPlist(self):
        'string representation: label and list of snps'
        
        snpString = ' '.join(snp.label for snp in self.snpList)
        return '%-25s %s' % (self.label, snpString)
    
    def strTreeTable(self):
        'string representation: line for tree table'
        
        label = self.label.split('/')[0]    # re-split to enable "Root"

        hgShortWithSNP = self.strHgShortWithSNP()
        if not hgShortWithSNP:              # for branches without SNPs
            hgShortWithSNP = '.'
        
        if self.isRoot():
            parentDFSrank, parentLabel = '.', '.'
        else:
            parentDFSrank = str(self.parent.DFSrank)
            parentLabel = self.parent.label.split('/')[0]
            
        return '%-5d %-25s %-20s %-6s %-25s' % \
            (self.DFSrank, label, hgShortWithSNP, parentDFSrank, parentLabel)
        
    @property
    def mostHighlyRankedSNP(self):
        'the most highly ranked SNP'
        
        return Node.mostHighlyRankedSNPonList(self.snpList)
    
    # static methods, including class variable setters
    #----------------------------------------------------------------------
    @staticmethod
    def setTreeConfigAndArgs(tree):
        'enables Node class to know about the tree instance, config, and args'
        
        Node.tree      = tree
        Node.config    = tree.config
        Node.args      = tree.args
        Node.errAndLog = tree.config.errAndLog
        if Node.args.writeContentMappings:
            Node.buildPageDict()

    @staticmethod
    def buildPageDict():
        'builds a dictionary of 23andMe content pages'
        
        with open(Node.config.pagesFN) as pagesFile:
            for line in pagesFile:
                haplogroup, hgSNP = line.strip().split()
                page = Page(haplogroup, hgSNP)
                Node.pageDict[haplogroup] = page
                if page.snpLabel:
                    Node.pageDict[page.snpLabel] = page

    @staticmethod
    def mostHighlyRankedSNPonList(snpList):
        '''returns the most highly ranked SNP on a list. 
            the purpose of this method is to centralize the knowledge that SNP lists 
            are now ranked in priority order rather than in reverse-priority order'''   
        
        if snpList:
            return snpList[0]
        else:
            return None
            
    # setters, mutaters
    #----------------------------------------------------------------------
    def setLabel(self, label):
        'sets label, haplogroup, and hgShort'
        
        self.label = label
        labelList = label.split('/')
        
        if self.isRoot():
            self.haplogroup = self.config.rootHaplogroup
            Node.tree.hg2nodeDict[self.haplogroup] = self
        else:
            self.haplogroup = labelList[0]
            self.setHgShort()
        
        for key in labelList:
            Node.tree.hg2nodeDict[key] = self
            
        if self.haplogroup in Node.pageDict:
            self.page = Node.pageDict[self.haplogroup]

    def setHgShort(self):
        'returns first 3-5 characters of specified haplogroups and first 2 of others'
        
        for numChars in xrange(5, 2, -1):
            if self.haplogroup[:numChars] in Node.config.mutiCharHgShortFormSet:
                self.hgShort = self.haplogroup[:numChars]
                break
        if not self.hgShort:
            self.hgShort = self.haplogroup[:2]

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
            page = Node.pageDict[snp.label]
            if self.page and self.page is not page:
                sys.exit('ERROR. Inconsistent pages assigned to node: %s\n' % self.label + \
                         'Page 1: %s\n' % self.page + \
                         'Page 2: %s\n' % page)
            self.page = page
        
    def sortSNPlist(self):
        'sorts snps by priority ranking'
    
        self.snpList = sorted(self.snpList, 
                              key=attrgetter('labelLettersRank', 
                                             'labelLetters', 
                                             'labelNumber'))
    
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
        
        genotypedSnpList = [snp for snp in self.snpList \
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
        for strLen in xrange(startLength, endLength):
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
        
        bfTreeFile.write('%s\n' % self)
        nodeDeque = deque(self.childList)
        while nodeDeque:
            node = nodeDeque.popleft()
            bfTreeFile.write('%s\n' % node)
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
    def writeNewick(self, newickFN, alignTips=False, platformVersion=None):
        'write Newick string for the subtree rooted at this node'
        
        treeString = '%s;\n' % self.buildNewickStringRecursive(alignTips, platformVersion)
        with open(newickFN, 'w') as outFile:
            outFile.write(treeString)
            
        if alignTips:
            treeDescriptor = 'aligned '
        elif platformVersion:
            treeDescriptor = 'platform v%d ' % platformVersion
        else:
            treeDescriptor = ''
        Node.errAndLog('Wrote %stree:\n    %s\n\n' % (treeDescriptor, newickFN))

    def buildNewickStringRecursive(self, alignTips=False, platformVersion=None):
        'recursively builds Newick string for the subtree rooted at this node'
        
        if not self.isLeaf():
            childStringList = list()
            for child in self.childList[::-1]:
                childString = child.buildNewickStringRecursive(alignTips, platformVersion)
                childStringList.append(childString)
            treeStringPart1 = '(%s)' % ','.join(childStringList)
        else:
            treeStringPart1 = ''

        branchLength = self.getBranchLength(alignTips, platformVersion)
        if alignTips:
            label = '%s:%d' % (self.label, branchLength)
        elif branchLength is None or (self.isLeaf() and branchLength == 0):
            label = self.label
        elif branchLength > 0:
            label = '%s|%d:%d' % (self.label, branchLength, branchLength)
        else:
            label = ':0.5'
        
        treeString = '%s%s' % (treeStringPart1, label)
        return treeString
