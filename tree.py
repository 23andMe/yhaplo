# David Poznik
# 2015.01.12
# tree.py
#
# Defines the Tree class.
#----------------------------------------------------------------------
from __future__ import absolute_import
import csv
import os
import re
import six
import sys
from collections import defaultdict, deque
from operator import attrgetter
from six.moves import range

from . import utils
from .node import Node
from .path import Path
from .snp import SNP, DroppedMarker

class Tree(object):
    '''
    A tree has single Node instance (the root) as well as a dictionary
    that maps haplogroup labels to node instances.
    '''
    
    def __init__(self, config):
        self.config = config
        self.args = config.args
        self.errAndLog = config.errAndLog
        self.maxDepth = 0
        
        # node info
        self.hg2nodeDict = dict()
        self.depthFirstNodeList = list()
        
        # snp info
        self.snpDict = dict()  # keys: name, (haplogroup, position), position
        self.snpList = list()
        self.snpPosSet = set()
        self.snpNameSet = set()
        self.preferredSNPnameSet = set()
        self.representativeSNPnameSet = set()
        self.multiAllelicOldPosSet = set()
        self.multiAllelicNewPosSet = set()
        self.isoggOmitSet = set()
        self.isoggCorrectionDict = dict()
        self.isoggCountsDict = defaultdict(int)
        self.numSNPsCorrected = 0
        
        # build tree
        self.root = self.buildTreeFromNewick()
        if self.args.primaryOnly:
            self.setDepthFirstNodeList()
        else:
            self.importIsoggSnps()
        self.setSearchRoot()
        self.optionalTraversalOutput()
        if self.args.writeContentMappings:
            self.writeContentMappings()
        
    # setters
    #----------------------------------------------------------------------
    def setSearchRoot(self):
        'set node from which to start haplogroup-calling traversals'
        
        if self.args.alternativeRoot:
            alternativeRootHg = self.args.alternativeRoot
            if alternativeRootHg in self.hg2nodeDict:
                self.searchRoot = self.hg2nodeDict[alternativeRootHg]
                self.errAndLog('Will start haplogroup assignment traversal from:\n' +
                               '    %s\n\n' % alternativeRootHg)
            else:
                sys.exit('\nERROR. Cannot start traversal ' +
                         'from non-existant haplogroup: %s\n' % alternativeRootHg)
        else:
            self.searchRoot = self.root
            
    def setDepthFirstNodeList(self):
        'build node list from depth-first pre-order traversal'
        
        self.depthFirstNodeList = self.root.getDepthFirstNodeList()
        for DFSrank, node in enumerate(self.depthFirstNodeList):
            node.setDFSrank(DFSrank)
        
    # traversals
    #----------------------------------------------------------------------
    def optionalTraversalOutput(self):
        'optional tree-traversal output'
        
        if self.args.traverseBF:
            self.writeBreadthFirst()
        if self.args.traverseDF:
            self.writeDepthFirstPreOrder()
        if self.args.writeTreeTable:
            self.writeTreeTable()
        if self.args.mrcaHaplogroupList:
            self.queryMRCA()
        if self.args.querySNPname:
            self.querySNPpath()

    def writeBreadthFirst(self):
        'writes bread-first traversal in pipe/dot format'
        
        bfTreeFN = (self.config.bfPrimaryTreeFN if self.args.primaryOnly
                    else self.config.bfTreeFN)
        with open(bfTreeFN, 'w') as bfTreeFile:
            self.root.writeBreadthFirstTraversal(bfTreeFile)
        
        self.errAndLog('Wrote breadth-first tree traveral:\n    %s\n\n' % bfTreeFN)
    
    def writeDepthFirstPreOrder(self):
        'writes depth-first pre-order traversal in pipe/dot format'
        
        dfTreeFN = (self.config.dfPrimaryTreeFN if self.args.primaryOnly
                    else self.config.dfTreeFN)
        with open(dfTreeFN, 'w') as dfTreeFile:
            for node in self.depthFirstNodeList:
                dfTreeFile.write('%s\n' % node.strDotPipeDepth())
                
        self.errAndLog('Wrote depth-first tree traveral:\n    %s\n\n' % dfTreeFN)
       
    def writeTreeTable(self):
        'writes depth-first pre-order traversal in table format'
        
        treeTableFN = self.config.treeTableFN
        headerList = '#index ycc_label label parent_index parent_label'.split()
        with open(treeTableFN, 'w') as treeTableFile:
            treeTableFile.write('%s\n' % '\t'.join(headerList))
            for node in self.depthFirstNodeList:
                treeTableFile.write('%s\n' % node.strTreeTableRow())
                
        self.errAndLog('Wrote tree table:\n    %s\n\n' % treeTableFN)
        
    def writeContentMappings(self):
        'writes the best 23andMe content page for each node'
        
        with open(self.config.pageUpdatesFN, 'w') as pageUpdatesFile:
            pageUpdatesFile.write('%-10s %-10s %-10s %s\n' %
                                  ('yccOld', 'SNP', 'hgSNP', 'ycc'))
            for page in Node.pageList:
                pageUpdatesFile.write('%s\n' % page.strFull())
            
        with open(self.config.pageMappingsFN, 'w') as pageMappingsFile:
            for node in self.depthFirstNodeList:
                ancestorNode = node
                while node.page is None:
                    ancestorNode = ancestorNode.parent
                    node.page = ancestorNode.page
                
                pageMappingsFile.write('%-25s %-15s | %s\n' %
                    (node.haplogroup, node.hgSNP if node.hgSNP else '.', node.page))
                
        self.errAndLog('%s23andMe content pages\n\n' % utils.DASHES +
            'Read  %4d titles:\n    %s\n\n' %
                (len(Node.pageList), self.config.pagesFN) +
            'Wrote %4d updates:\n    %s\n\n' %
                (len(Node.pageList), self.config.pageUpdatesFN) +
            'Wrote %4d mappings:\n    %s\n\n' %
                (len(self.depthFirstNodeList), self.config.pageMappingsFN))

    def queryMRCA(self):
        'writes MRCA of two haplogroups'
        
        mrcaHaplogroupList = self.args.mrcaHaplogroupList
        
        if type(mrcaHaplogroupList) != list or len(mrcaHaplogroupList) != 2:
            sys.exit('ERROR. mrca expects a list of 2 haplogroups: %s\n' %
                     mrcaHaplogroupList)
        haplogroup1, haplogroup2 = mrcaHaplogroupList
        node1 = self.hg2nodeDict[haplogroup1]
        node2 = self.hg2nodeDict[haplogroup2]
        mrca = node1.mrca(node2)
        self.errAndLog('%sMRCA Query\n\n' % utils.DASHES +
                       'Haplogroup 1: %s\n' % node1.haplogroup +
                       'Haplogroup 2: %s\n' % node2.haplogroup +
                       'MRCA: %s\n\n' % mrca.haplogroup)
        
    def querySNPpath(self):
        'lists phylogenetic path for a query SNP'
        
        queryName = self.args.querySNPname
        self.errAndLog('%sSNP Query: %s\n\n' % (utils.DASHES, queryName))
        snp = self.snpDict.get(queryName, None)
        
        if snp:
            for node in snp.backTracePath():
                self.errAndLog('%s\n' % node.strSimple())
            if snp.label != queryName:
                self.errAndLog('\nNote: %s is an alias of %s.\n' % (queryName, snp.label))
        else:
            self.errAndLog('Not found.\n')

        self.errAndLog('\n')
        
    # write newick files
    #----------------------------------------------------------------------
    def writeNewick(self):
        'write tree as is and with aligned terminal branch lengths'

        if not self.config.suppressOutputAndLog:
            self.errAndLog('%sWriting trees...\n\n' % utils.DASHES)
            self.root.writeNewick(self.config.yccTreeFN)
            self.root.writeNewick(self.config.hgsnpTreeFN, useHgSNPlabel=True)
            self.root.writeNewick(self.config.alignedYccTreeFN, alignTips=True)
            self.root.writeNewick(self.config.alignedHgsnpTreeFN,
                                  useHgSNPlabel=True, alignTips=True)
            if self.args.writePlatformTrees:
                self.writePlatformTrees()
    
    def writePlatformTrees(self):
        'write trees whose branch lengths are numbers of platform sites'
        
        for platformVersion in range(1, self.config.maxPlatformVersionPlusOne):
            self.root.writeNewick(self.config.platformYccTreeFNtp % platformVersion,
                                  platformVersion=platformVersion)
            self.root.writeNewick(self.config.platformHgsnpTreeFNtp % platformVersion,
                                  useHgSNPlabel=True, platformVersion=platformVersion)

    # query
    #----------------------------------------------------------------------
    def identifyPhylogeneticPath(self, sample):
        '''
        conducts a modified breadth-first search (bfs) to identify
        the phylogenetic path leading from the root to the most terminal branch
        representing a sample's haplogroup.
        
        returns: best path, list of SNPs observed in the ancestral state.

        key differences from a standard bfs are:
        - stopping condition is robust to genotype error, homoplasy, etc.
        - collapsing condition to speed up and (marginally) improve accuracy
        
        when the stopping condition is met, adds the current path to a list.
        at the end, post-processes this list and selects the best element.
        
        the stopping condition is a disjunction of three literals.
        the first is trivial:
        
        a. node.isLeaf()
            we cannot go any further

        the following table enumerates possible cases for the two other literals.
        #anc: number of ancestral alleles observed on a branch
        #der: number of derived alleles observed on the branch
                  (only considered if #anc == 2)
        stop: whether or not to stop

        | #anc | #der | stop | reason
        |------|------|------|--------------------------------------------------------
        | 0, 1 |    . |   no | insufficient evidence to stop, or none at all
        |    2 |   1+ |   no | given evidence to continue, do so for robustness
        |    2 |    0 |  yes | reasonable evidence to stop and no evidence to continue
        |   3+ |    . |  yes | strong evidence to stop

        b. row 4: numAncestral > 2
           der =  0: compelling evidence to stop.
           der = 1+: the sample's lineage probably diverges from the known tree here.

        c. row 3: numAncestral == 2 and numDerived == 0
            it is safe to assume that this path will not yield fruit.

        these conditions are robust to the most challenging case:
        when just a single SNP is genotyped on a branch, and the observed genotype
        corresponds to the ancestral allele due to genotype error, homoplasy,
        or an undetected isogg error. when at least one derived allele is observed,
        the conditions are also robust to two false ancestral alleles on a branch.
        '''
        
        pathDeque = Path.createPathDeque(self.searchRoot.childList)
        stoppedPathList = list()
        ancSNPfullList = list()
        while pathDeque:
            path = pathDeque.popleft()
            ancSNPlist, derSNPlist = path.node.assessGenotypes(sample)
            path.updateWithBranchAssessment(ancSNPlist, derSNPlist)
            ancSNPfullList.extend(ancSNPlist)
            numAncestral, numDerived = len(ancSNPlist), len(derSNPlist)

            if self.args.writeAncDerCounts:
                sample.appendAncDerCountTuple(path.node, numAncestral, numDerived)
            
            if (path.node.isLeaf()
                or (numAncestral  > self.config.args.ancStopThresh)
                or (numAncestral == self.config.args.ancStopThresh and numDerived == 0)):
                stoppedPathList.append(path)
            else:
                if numDerived >= self.config.args.derCollapseThresh:
                    pathDeque = deque()
                pathDeque.extend(path.fork(path.node.childList))
        
        bestPath = Path.postProcessPathListAndSelectBest(stoppedPathList)
        return bestPath, ancSNPfullList
    
    def getDFSrank(self, haplogroup):
        'returns the DFS rank of a haplogroup'
        
        return self.hg2nodeDict[haplogroup].DFSrank
    
    # build tree from Newick-formatted text file
    #----------------------------------------------------------------------
    def buildTreeFromNewick(self):
        '''
        Reads in a Newick-formatted tree, strips out bootstraps,
        tokenizes it, and initiates tree building.
        Returns a node instance: the root.
        '''
        
        utils.checkFileExistence(self.config.primaryTreeFN, 'Primary tree')
        with open(self.config.primaryTreeFN, 'r') as treeFile:
            treeString = treeFile.readline().strip()
        self.errAndLog('\n%sRead primary tree topology:\n    %s\n\n' %
                       (utils.DASHES, self.config.primaryTreeFN))

        '''
        Tokenization:
            a. strip out bootstraps: text within brackets
            b. split on any semantic token: [%s]
            c. but group to retain retain tokens themselves: ()
            d. then drop empty tokens from splitting adjacent semantic tokens
        '''
        treeString = re.subn(r'\[.*?\]', '', treeString)[0]
        treeList = re.split('([%s])' % self.config.newickSemanticTokenString, treeString)
        treeList = [token for token in treeList if token is not '']
        treeDeque = deque(treeList)
    
        hasLengths = ':' in treeDeque  # determine whether tree has lengths
        root = self.addChildSubtreeFromNewickDeque(None, treeDeque, hasLengths)
        root.writeNewick(self.config.alignedPrimaryTreeFN, alignTips=True)
        
        return root
    
    def addChildSubtreeFromNewickDeque(self, parent, treeDeque, hasLengths):
        '''
        Recursively processes a deque of Newick tokens to build a tree.
        Each call constructs one subtree and returns its root.
        1. Recursive case: an open paren indicates a compound subtree.
            The function calls itself to add the first child.
        2. Base case: an alphanumeric label indicates a leaf.
            Return a simple leaf node.
        3. Following the first child subtree, there will be
            an arbitrary number of sibling subtrees, each preceeded by a comma.
            The function calls itself to add each in turn.
        4. The end of a subtree signaled by a close paren.
            At this point, add a label and/or length, if either are provided.
        '''
        
        #-------------------------------------------------------------------------
        # first node of subtree
        node = Node(parent = parent, tree = self)
        token = treeDeque.popleft()
        if token == '(':    # recursive case: compound subtree
            self.addChildSubtreeFromNewickDeque(node, treeDeque, hasLengths)
        else:               # base case: leaf tree
            node.setLabel(token)
            if hasLengths:
                Tree.processNewickLength(node, treeDeque)
            return node
    
        #-------------------------------------------------------------------------
        # second through nth nodes of subtree
        token = treeDeque.popleft()
        while token == ',':
            self.addChildSubtreeFromNewickDeque(node, treeDeque, hasLengths)
            token = treeDeque.popleft()
    
        #-------------------------------------------------------------------------
        # end of subtree
        Tree.verifyToken(token, ')')
        node.reverseChildren()
        token = treeDeque.popleft()
        if token not in self.config.newickSemanticTokenSet:
            node.setLabel(token)
        if hasLengths and treeDeque[0] != ';':
            self.config.processNewickLength(node, treeDeque)
        return node

    @staticmethod
    def verifyToken(observed, expected):
        'exits program if observed and expected strings do not match'
        
        if observed != expected:
            sys.exit('ERROR. Malformed newick file.\n' +
                     'Expected this token: %s\n' % expected +
                     'Got this one:        %s\n' % observed)
    
    @staticmethod
    def processNewickLength(node, treeDeque):
        'processes a Newick-format branch length of the form :length'
        
        Tree.verifyToken(treeDeque.popleft(), ':')    # next token should be colon
        branchLength = float(treeDeque.popleft())     # then branch length
        node.setBranchLength(branchLength)

    # import SNPs and assign to branches
    #----------------------------------------------------------------------
    def importIsoggSnps(self):
        'import ISOGG SNPs'
        
        SNP.setClassVariables(self)
        self.readPreferredSNPnameSet()
        self.readRepresentativeSNPnameSet()
        self.readIsoggMultiAllelicPosSet()
        self.readIsoggOmitSet()
        self.readIsoggCorrectionsDict()
        self.parseIsoggTable()
        self.setDepthFirstNodeList()
        self.sortSNPlistsAndSetRepresentatives()
        self.writeIsoggCounts()
        self.writeUniqueSNPtable()
        self.writeNewick()
        self.checkMultiAllelics()
    
    def readPreferredSNPnameSet(self):
        '''reads a set of widely known SNP names. presence on this list is
            the primary selection criterion for SNP labels'''
        
        preferredSNPnamesFN = self.config.preferredSNPnamesFN
        
        utils.checkFileExistence(preferredSNPnamesFN, 'Preferred SNP names')
        with open(preferredSNPnamesFN, 'r') as preferredSNPnamesFile:
            for line in preferredSNPnamesFile:
                self.preferredSNPnameSet.add(line.strip())
                
        self.errAndLog(
            '%sRead preferred SNP names\n' % utils.DASHES +
            '%6d SNP names: %s\n\n' %
                (len(self.preferredSNPnameSet), preferredSNPnamesFN))
    
    def readRepresentativeSNPnameSet(self):
        'reads the names of SNPs deemed representative for their respective lineages'
        
        isoggRepSNPfn = self.config.isoggRepSNPfn
        otherRepSNPfn = self.config.otherRepSNPfn
        countsDicts = defaultdict(int)
        
        set1 = set()
        utils.checkFileExistence(isoggRepSNPfn, 'First representative SNPs')
        with open(isoggRepSNPfn, 'r') as isoggRepSNPfile:
            for line in isoggRepSNPfile:
                countsDicts['lines'] += 1
                snpAliasesString = line.strip().split()[1]
                if snpAliasesString != '.':
                    countsDicts['haplogroups'] += 1
                    for snpAliases in snpAliasesString.split(','):
                        countsDicts['snps'] += 1
                        for snpName in snpAliases.split('/'):
                            set1.add(snpName)
        
        set2 = set()
        utils.checkFileExistence(otherRepSNPfn, 'Second representative SNPs')
        with open(otherRepSNPfn, 'r') as otherRepSNPfile:
            for line in otherRepSNPfile:
                set2.add(line.strip().split()[1])
        
        self.representativeSNPnameSet = set1 | set2
        self.errAndLog(
            'Read representative SNPs\n' +
            '%6d haplogroups in: %s\n' % (countsDicts['lines'], isoggRepSNPfn) +
            '%6d haplogroups with at least one ISOGG-designated representative SNP\n' %
                countsDicts['haplogroups'] +
            '%6d SNPs, as some haplogroups have more than one representative\n' %
                countsDicts['snps'] +
            '%6d SNP names, including aliases\n' % len(set1) +
            '%6d additional representative SNPs read from: %s\n' % (len(set2), otherRepSNPfn) +
            '%6d total SNP names\n\n' % len(self.representativeSNPnameSet))
    
    def readIsoggMultiAllelicPosSet(self):
        'reads list of positions to exclude because of multiple alleles'
        
        if not os.path.isfile(self.config.isoggMultiAllelicFN):
            return
        with open(self.config.isoggMultiAllelicFN, 'r') as multiFile:
            for line in multiFile:
                position = int(line.strip())
                self.multiAllelicOldPosSet.add(position)
        
    def readIsoggOmitSet(self):
        'reads a list of SNPs to omit from ISOGG db'
        
        for isoggOmitFN in self.config.isoggOmitFNlist:
            if not os.path.isfile(isoggOmitFN):
                continue
            with open(isoggOmitFN, 'r') as omitFile:
                for line in omitFile:
                    lineList = line.strip().split()
                    if len(lineList) > 0 and lineList[0] != '#':
                        position, mutation = lineList[2:4]
                        self.isoggOmitSet.add((position, mutation))
        
    def readIsoggCorrectionsDict(self):
        'reads a list of SNPs to correct from ISOGG db'
        
        for isoggCorrectionsFN in self.config.isoggCorrectionsFNlist:
            if not os.path.isfile(isoggCorrectionsFN):
                continue
            with open(isoggCorrectionsFN, 'r') as correctionsFile:
                for line in correctionsFile:
                    lineList = line.strip().split()
                    if len(lineList) > 0 and lineList[0] != '#':
                        haplogroup, position, mutation, aliases = lineList[1:5]
                        for alias in aliases.split(','):
                            self.isoggCorrectionDict[alias] = (
                                (haplogroup, position, mutation))
    
    def parseIsoggTable(self):
        'parses ISOGG table'
        
        # input reader
        utils.checkFileExistence(self.config.isoggFN, 'Isogg')
        isoggInFile = open(self.config.isoggFN, 'r')
        isoggReader = csv.reader(isoggInFile, delimiter='\t')
        next(isoggReader) # ignore header
        
        # output file handles
        if self.config.suppressOutputAndLog:
            isoggOutFile = None
            isoggDropOutFile = None
        else:
            isoggOutFile = open(self.config.cleanedIsoggFN, 'w')
            isoggDropOutFile = open(self.config.droppedIsoggFN, 'w')
        
        droppedMarkerList = list()
        
        for lineList in isoggReader:
            self.isoggCountsDict['read'] += 1

            # clean up data row and extract values
            lineList = [element.strip() for element in lineList]
            if lineList[1] == '':   # when present, remove extra tab after snp name
                del lineList[1]
            if len(lineList) != 6:
                self.isoggCountsDict['badLines'] += 1
                continue
            name, haplogroup, _, _, position, mutation = lineList
            
            # apply corrections
            if name in self.isoggCorrectionDict:
                haplogroup, position, mutation = self.isoggCorrectionDict[name]
                self.numSNPsCorrected += 1
                
            # identify markers to drop
            recordIsBad, markerIsOkToRepresentNode = (
                self.checkIsoggRecord(name, haplogroup, position, mutation))
            if recordIsBad:
                self.isoggCountsDict['dropped'] += 1
                if isoggDropOutFile:
                    isogg_drop_output = six.ensure_str('%-10s %-25s %8s %s\n' %
                        (six.ensure_text(name), haplogroup, position, mutation))
                    isoggDropOutFile.write(isogg_drop_output)
                if markerIsOkToRepresentNode:
                    droppedMarkerList.append(DroppedMarker(name, haplogroup))
                continue
            
            # process retained SNPs
            self.isoggCountsDict['retained'] += 1
            position = int(position)
            if isoggOutFile:
                isoggOutFile.write('%-10s %-25s %8d %s\n' %
                                   (name, haplogroup, position, mutation))
            self.constructSNP(name, haplogroup, position, mutation)
        
        self.addDroppedMarkersToNodes(droppedMarkerList)
        utils.closeFiles([isoggInFile, isoggOutFile, isoggDropOutFile])

    def constructSNP(self, name, haplogroup, position, mutation):
        '''
        typically, instantiates a SNP and adds it to various containers.
        note that when SNPs are instantiated, they are added to the tree,
        and this process may entail growing the tree to include the corresponding node.
        
        more specialized things occur if a SNP already exists at this position.
        '''
        
        if self.hg2nodeDict:
            ancestral, derived = mutation[0], mutation[3]
            snpKey = (haplogroup, position)
            
            if snpKey in self.snpDict:      # snp exists under an alias
                snp = self.snpDict[snpKey]
                if snp.isAncestral(ancestral) and snp.isDerived(derived):
                    snp.addName(name)
                    self.snpDict[name] = snp
                else:
                    newSNP = SNP(name, haplogroup, position, ancestral, derived)
                    sys.exit('\n\nERROR! Conlicting SNPs:\n%s\n%s\n' %
                             (snp, newSNP))
            else:
                if position in self.snpDict: # another snp with same position
                    oldSNP = self.snpDict[position]
                    if (ancestral not in oldSNP.alleleSet
                        or derived not in oldSNP.alleleSet):
                        self.multiAllelicNewPosSet.add(position)
                
                # typical behavior
                snp = SNP(name, haplogroup, position, ancestral, derived)
                self.snpDict[(haplogroup, position)] = snp
                self.snpDict[name] = snp
                self.snpDict[position] = snp
                self.snpList.append(snp)
                self.snpPosSet.add(position)
                self.snpNameSet.add(name)
                self.isoggCountsDict['unique'] += 1
                
    def addDroppedMarkersToNodes(self, droppedMarkerList):
        'adds dropped markers to coresponding nodes'
        
        for droppedMarker in droppedMarkerList:
            droppedMarker.addToNode()
        
    def sortSNPlistsAndSetRepresentatives(self):
        'for each node, sorts snps by priority ranking and selects the best representative'
        
        if not self.depthFirstNodeList:
            self.setDepthFirstNodeList
            
        for node in self.depthFirstNodeList:
            node.prioritySortSNPlistAndSetHgSNP()
        
    def writeIsoggCounts(self):
        config = self.config
        countsDict = self.isoggCountsDict
        numAltNames = countsDict['retained'] - countsDict['unique']
        
        logTextList = [
            '%sRead ISOGG SNP data:\n    %s\n' % (utils.DASHES, config.isoggFN),
            '  %5d SNPs read'                          % countsDict['read'],
            '  %5d corrected '                         % self.numSNPsCorrected +
            'based on:\n            %s\n' % ('\n'+' '*12).join(config.isoggCorrectionsFNlist),
            '- %5d SNPs dropped'                       % countsDict['dropped'],
            '        %5d flagged as not meeting quality guidelines' % countsDict['qc'],
            '        %5d tree location approximate'    % countsDict['approxLoc'],
            '        %5d removed, flagged as provisional, or otherwise problematic' %
                                                         countsDict['provisional'],
            '        %5d non-SNPs'                     % countsDict['nonSNP'],
            '        %5d excluded as multiallelic '    % countsDict['multiallelic'] +
            'based on:\n                  %s'          % config.isoggMultiAllelicFN,
            '        %5d duplicated names'             % countsDict['duplicatedNames'],
            '        %5d explicitly excluded '         % countsDict['omitted'] +
            'based on:\n                  %s' % ('\n'+' '*18).join(config.isoggOmitFNlist),
            '- %5d bad lines'                          % countsDict['badLines'],
            '= %5d SNPs retained\n'                    % countsDict['retained'],
            '- %5d alternative names'                  % numAltNames,
            '= %5d unique SNPs added to the tree\n'    % countsDict['unique'],
        ]
        
        if not config.suppressOutputAndLog:
            logTextList.extend([
                'Wrote summary tables',
                '    dropped:  %s'   % config.droppedIsoggFN,
                '    retained: %s'   % config.cleanedIsoggFN,
                '    unique:   %s\n' % config.uniqueIsoggFN,
            ])
        
        self.errAndLog(('\n' + ' '*4).join(logTextList) + '\n')

    def writeUniqueSNPtable(self):
        'sort unique snp list by phylogeny and position; write to file'
        
        if not self.config.suppressOutputAndLog:
            self.snpList = sorted(self.snpList, key = attrgetter('DFSrank', 'position'))
            with open(self.config.uniqueIsoggFN, 'w') as uniqueIsoggFile:
                for snp in self.snpList:
                    uniqueIsoggFile.write('%s\n' % snp.strWithAllNames())
    
    def checkIsoggRecord(self, name, haplogroup, position, mutation):
        '''
        returns a tuple of booleans: (recordIsBad, nameIsOkToRepresentNode)
        
        how to interpret TRUE values:
        1. recordIsBad               : do not use this marker for classification
        2. markerIsOkToRepresentNode : if no SNPs are retained for the corresponding node,
                it's OK to use this marker name for the node's hgSNP representation
        '''
        
        if name.endswith('^'):
            self.isoggCountsDict['qc'] += 1
            return True, True
        
        if haplogroup.find('~') >= 0:
            self.isoggCountsDict['approxLoc'] += 1
            return True, False      # second value irrelevant: no corresponding node
    
        if (haplogroup.find('Investigation') >= 0
                or haplogroup.find('Notes') >= 0
                or haplogroup.find('Private') >= 0
                or haplogroup.find('Removed') >= 0
                or haplogroup.find('Withdrawn') >= 0
                or haplogroup.find('Freq. Mut.') >= 0
                or len(haplogroup) < 1):
            self.isoggCountsDict['provisional'] += 1
            return True, False      # second value irrelevant: no corresponding node
    
        if (len(mutation) != 4
                or mutation.find('?') >= 0
                or position.find('..') >= 0):
            self.isoggCountsDict['nonSNP'] += 1
            return True, True
        
        if int(position) in self.multiAllelicOldPosSet:
            self.isoggCountsDict['multiallelic'] += 1
            return True, True
        
        if name in self.snpNameSet:
            self.isoggCountsDict['duplicatedNames'] += 1
            return True, False
        
        if (position, mutation) in self.isoggOmitSet:
            self.isoggCountsDict['omitted'] += 1
            return True, False
    
        try:
            position = int(position)
        except ValueError:
            self.errAndLog('\nERROR. Invalid position: %s\n' % position)
            return True, False
        
        return False, True

    def checkMultiAllelics(self):
        'checks for mutliallelic variants and writes list to file'
        
        if len(self.multiAllelicNewPosSet) > 0:
            if not self.config.suppressOutputAndLog:
                with open(self.config.multiAllelicFoundFN, 'w') as outFile:
                    for position in sorted(list(self.multiAllelicNewPosSet)):
                        outFile.write('%8d\n'% position)
                    
            self.errAndLog('%s*** Dectected %d multiallelic positions. ***\n\n' %
                                (utils.DASHES, len(self.multiAllelicNewPosSet)) +
                           'Please do the following and then re-run:\n' +
                           '    cat %s >> %s\n\n' %
                                (self.config.multiAllelicFoundFN,
                                 self.config.isoggMultiAllelicFN))

    def findOrCreateNode(self, haplogroup):
        '''
        given a haplogroup, returns corresponding node if it exists.
        if not, first serially split most recent ancestor that does exists
        until there is a place for it.
        
        Note: the loop has to start at 1, because myString[:0] --> ''
        '''
        
        if haplogroup in self.hg2nodeDict:
            return self.hg2nodeDict[haplogroup]
        
        for numCharsToChop in range(1, len(haplogroup)):
            ancestorString = haplogroup[:-numCharsToChop]
            if ancestorString in self.hg2nodeDict:
                ancestor = self.hg2nodeDict[ancestorString]
                return ancestor.serialSplit(haplogroup)
        
        self.errAndLog('Unplaceable haplogroup: %s' % haplogroup)
