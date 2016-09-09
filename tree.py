# David Poznik
# 2015.01.12
# tree.py
#
# Defines the Tree class.
#----------------------------------------------------------------------
import csv
import os
import re
import sys
from collections import defaultdict, deque
from operator import methodcaller, attrgetter

import utils
from node import Node
from path import Path
from snp import SNP

class Tree(object):
    '''
    A tree has single Node instance (the root) as well as a dictionary
    that maps haplogroup labels to node instances.
    '''
    
    def __init__(self, config):
        self.config    = config
        self.args      = config.args
        self.errAndLog = config.errAndLog
        self.maxDepth  = 0
        
        # node info
        self.hg2nodeDict = dict()
        self.depthFirstNodeList = list()
        
        # snp info
        self.snpDict = dict()  # keys: name, (haplogroup, position), position 
        self.snpList               = list()
        self.snpPosSet             = set()
        self.snpNameSet            = set()
        self.multiAllelicOldPosSet = set()
        self.multiAllelicNewPosSet = set()
        self.isoggOmitSet          = set()
        self.isoggCorrectionDict   = dict()
        self.isoggCountsDict       = defaultdict(int)
        self.numSNPsCorrected      = 0
        
        # build tree
        self.root = self.buildTreeFromNewick()
        if self.args.primaryOnly:
            self.setDepthFirstNodeList()
        else:
            self.importIsoggSnps()
        self.setSearchRoot()
        self.optionalTraversalOutput()
        if Node.pageDict:
            self.setContentPages()
        
    # setters
    #----------------------------------------------------------------------
    def setSearchRoot(self):
        'set node from which to start haplogroup-calling traversals'
        
        if self.args.alternativeRoot:
            alternativeRootHg = self.args.alternativeRoot
            if alternativeRootHg in self.hg2nodeDict:
                self.searchRoot = self.hg2nodeDict[alternativeRootHg]
                self.errAndLog('Will start haplogroup assignment traversal from:\n' + \
                               '    %s\n\n' % alternativeRootHg)
            else:
                sys.exit('\nERROR. Cannot start traversal ' + \
                         'from non-existant haplogroup: %s\n' % alternativeRootHg)
        else:
            self.searchRoot = self.root
            
    def setDepthFirstNodeList(self):
        'build node list from depth-first pre-order traversal'
        
        self.depthFirstNodeList = self.root.getDepthFirstNodeList()
        for DFSrank, node in enumerate(self.depthFirstNodeList):
            node.setDFSrank(DFSrank)
            
    def setContentPages(self):
        'sets the 23andMe content page for each node'
        
        with open(self.config.pageMappingsFN, 'w') as pageMappingsFile:
            # coded to work for arbitrary traversals, but since this is depth-first,
            # each page will be assigned in at most one step 
            
            for node in self.depthFirstNodeList:
                ancestorNode = node
                while node.page is None:
                    ancestorNode = ancestorNode.parent
                    node.page = ancestorNode.page
                
                pageMappingsFile.write('%-25s %s\n' % (node.haplogroup, node.page))
                
        self.errAndLog('%sWrote %d 23andMe content page mappings:\n    %s\n\n' % \
                       (utils.DASHES, len(self.depthFirstNodeList), 
                        self.config.pageMappingsFN))
        
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
            self.printMRCA(self.args.mrcaHaplogroupList)

    def writeBreadthFirst(self):
        'writes bread-first traversal in pipe/dot format'
        
        bfTreeFN = self.config.bfPrimaryTreeFN if self.args.primaryOnly \
              else self.config.bfTreeFN
        with open(bfTreeFN, 'w') as bfTreeFile:
            self.root.writeBreadthFirstTraversal(bfTreeFile)
        
        self.errAndLog('Wrote breadth-first tree traveral:\n    %s\n\n' % bfTreeFN)
    
    def writeDepthFirstPreOrder(self):
        'writes depth-first pre-order traversal in pipe/dot format'
        
        dfTreeFN = self.config.dfPrimaryTreeFN if self.args.primaryOnly \
              else self.config.dfTreeFN
        with open(dfTreeFN, 'w') as dfTreeFile:
            for node in self.depthFirstNodeList:
                dfTreeFile.write('%s\n' % node)
                
        self.errAndLog('Wrote depth-first tree traveral:\n    %s\n\n' % dfTreeFN)
       
    def writeTreeTable(self):
        'writes depth-first pre-order traversal in table format'
        
        treeTableFN = self.config.treeTableFN
        with open(treeTableFN, 'w') as treeTableFile:
            treeTableFile.write('%-5s %-25s %-20s %-6s %-25s\n' % \
                ('Index', 'Label', 'HgShort', 'Parent', 'ParentLabel'))
            for node in self.depthFirstNodeList:
                treeTableFile.write('%s\n' % node.strTreeTable())
                
        self.errAndLog('Wrote tree table:\n    %s\n\n' % treeTableFN)
        
    def printMRCA(self, mrcaHaplogroupList):
        'prints to stdout MRCA of two haplogroups'
        
        if type(mrcaHaplogroupList) != list or len(mrcaHaplogroupList) != 2:
            sys.exit('ERROR. mrca expects a list of 2 haplogroups: %s\n' % \
                     mrcaHaplogroupList)
        haplogroup1, haplogroup2 = mrcaHaplogroupList
        node1 = self.hg2nodeDict[haplogroup1]
        node2 = self.hg2nodeDict[haplogroup2]
        mrca = node1.mrca(node2)
        self.errAndLog('%sMRCA Query\n\n' % (utils.DASHES) + \
                       'Haplogroup 1: %s\n' % (node1.haplogroup) + \
                       'Haplogroup 2: %s\n' % (node2.haplogroup) + \
                       'MRCA: %s\n\n' % (mrca.haplogroup))
        
    # write newick files
    #----------------------------------------------------------------------
    def writeNewick(self):
        'write tree as is and with aligned terminal branch lengths'

        self.errAndLog('%sWriting trees...\n\n' % utils.DASHES)
        self.root.writeNewick(self.config.treeFN)
        self.root.writeNewick(self.config.alignedTreeFN, alignTips=True)
        if self.args.writePlatformTrees:
            self.writePlatformTrees()
    
    def writePlatformTrees(self):
        'write trees whose branch lengths are numbers of platform sites'
        
        for platformVersion in xrange(1, self.config.maxPlatformVersionPlusOne):
            newickFN = self.config.platformTreeFNtp % platformVersion
            self.root.writeNewick(newickFN, platformVersion = platformVersion)

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
            
            if path.node.isLeaf() \
                or (numAncestral  > self.config.ancStopThresh) \
                or (numAncestral == self.config.ancStopThresh and numDerived == 0):
                stoppedPathList.append(path)
            else:
                if numDerived >= self.config.derCollapseThresh:
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
        
        primaryTreeFN = self.config.primaryTreeFN
        if not os.path.isfile(primaryTreeFN):
            sys.exit('ERROR. Primary tree file not found: %s\n' % primaryTreeFN)
        with open(primaryTreeFN, 'r') as treeFile:
            treeString = treeFile.readline().strip()
        self.errAndLog('\n%sRead primary tree topology:\n    %s\n\n' % \
                       (utils.DASHES, primaryTreeFN))

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
            sys.exit('''
            ERROR. Malformed newick file.
            Expected this token: %s
            Got this one:        %s
            ''' % (expected, observed))
    
    @staticmethod
    def processNewickLength(node, treeDeque):
        'processes a Newick-format branch length of the form :length'
        
        Tree.verifyToken(treeDeque.popleft(), ':')    # next token should be colon
        branchLength = float(treeDeque.popleft())     # then branch length
        node.setBranchLength(branchLength)

    # import SNPs and assign to branches
    #----------------------------------------------------------------------
    def importIsoggSnps(self):
        '''1. set up SNP class variables.
           2. read ISOGG corrections and omissions.
           3. set up file handles for table parsing.
           4. call ISOGG table parser.
           5. close files.
           6. post-processing.'''
        
        SNP.setClassVariables(self)
        
        self.readIsoggMultiAllelicPosSet()
        self.readIsoggOmitSet()
        self.readIsoggCorrectionsDict()
        
        isoggInFile      = open(self.config.isoggFN, 'r')
        isoggReader      = csv.reader(isoggInFile, delimiter='\t')
        isoggOutFile     = open(self.config.cleanedIsoggFN, 'w')
        isoggDropOutFile = open(self.config.droppedIsoggFN, 'w')
        
        self.parseIsoggTable(isoggReader, isoggOutFile, isoggDropOutFile)
    
        for File in [isoggInFile, isoggOutFile, isoggDropOutFile]:
            File.close()
            
        self.setDepthFirstNodeList()
        self.sortSNPlists()
        self.writeIsoggCounts()
        self.writeUniqueSNPtable()
        self.writeNewick()
        self.checkMultiAllelics()
    
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
                            self.isoggCorrectionDict[alias] = \
                                (haplogroup, position, mutation)
    
    def parseIsoggTable(self, isoggReader, isoggOutFile, isoggDropOutFile):
        'parses ISOGG table and adds SNPs to the tree'
        
        isoggReader.next()              # ignore header
        for lineList in isoggReader:
            self.isoggCountsDict['read'] += 1
            
            lineList = [element.strip() for element in lineList]
            if lineList[1] == '':   # extra tab after snp name
                del lineList[1]
            if len(lineList) != 6:
                self.isoggCountsDict['badLines'] += 1
                continue
            name, haplogroup, _, _, position, mutation = lineList
            if name in self.isoggCorrectionDict:
                haplogroup, position, mutation = self.isoggCorrectionDict[name]
                self.numSNPsCorrected += 1
            if self.isoggRecordIsBad(name, haplogroup, position, mutation):
                self.isoggCountsDict['dropped'] += 1
                isoggDropOutFile.write('%-10s %-25s %8s %s\n' % \
                    (name, haplogroup, position, mutation))
                continue
            position = int(position)
            self.isoggCountsDict['retained'] += 1
            isoggOutFile.write('%-10s %-25s %8d %s\n' % \
                (name, haplogroup, position, mutation))
            
            if len(self.hg2nodeDict) > 0:
                ancestral, derived = mutation[0], mutation[3]
                snpKey = (haplogroup, position)
                
                if snpKey in self.snpDict:      # snp exists under an alias
                    snp = self.snpDict[snpKey]
                    if snp.isAncestral(ancestral) and snp.isDerived(derived):
                        snp.addName(name)
                        self.snpDict[name] = snp
                    else:
                        newSNP = SNP(name, haplogroup, position, ancestral, derived)
                        sys.exit('\n\nERROR! Conlicting SNPs:\n%s\n%s\n' % \
                                 (snp, newSNP))
                else:
                    if position in self.snpDict: # another snp with same position
                        oldSNP = self.snpDict[position]
                        if ancestral not in oldSNP.alleleSet or \
                           derived not in oldSNP.alleleSet:
                            self.multiAllelicNewPosSet.add(position)
                            
                    snp = SNP(name, haplogroup, position, ancestral, derived)
                    self.snpDict[(haplogroup, position)] = snp
                    self.snpDict[name] = snp
                    self.snpDict[position] = snp
                    self.snpList.append(snp)
                    self.snpPosSet.add(position)
                    self.snpNameSet.add(name)
                    self.isoggCountsDict['unique'] += 1
        
    def sortSNPlists(self):
        '''sorts snps in reverse priority ranking so that the most derived SNP 
            observed for each sample best tags his haplogroup'''
        
        if self.depthFirstNodeList:
            for node in self.depthFirstNodeList:
                node.sortSNPlist()
        else:
            sys.exit('ERROR. Must build depthFirstNodeList ' + \
                     'before calling sortSNPlists.\n\n')
        
    def writeIsoggCounts(self):
        config     = self.config
        countsDict = self.isoggCountsDict
        numAltNames = countsDict['retained'] - countsDict['unique']
        
        self.errAndLog('''%sRead ISOGG SNP data:\n    %s\n
        %5d SNPs read
        %5d corrected based on:\n              %s\n
        %5d SNPs dropped and written:\n              %s
              %5d flagged as not meeting quality guidelines
              %5d tree location approximate
              %5d removed, flagged as provisional, or otherwise problematic
              %5d non-SNPs
              %5d excluded due to multiallelic positions based on:\n                    %s
              %5d duplicated names
              %5d explicitly excluded based on:\n                    %s
        %5d bad lines\n
        %5d SNPs retained and written:\n              %s
        %5d alternative names\n
        %5d unique SNPs added to the tree and written:\n              %s\n\n''' \
            % (utils.DASHES, config.isoggFN, 
                countsDict['read'], 
                self.numSNPsCorrected, ('\n'+' '*14).join(config.isoggCorrectionsFNlist),
                countsDict['dropped'], config.droppedIsoggFN,
                    countsDict['qc'], 
                    countsDict['approxLoc'], 
                    countsDict['provisional'], 
                    countsDict['nonSNP'], 
                    countsDict['multiallelic'], config.isoggMultiAllelicFN,
                    countsDict['duplicatedNames'], 
                    countsDict['omitted'], ('\n'+' '*20).join(config.isoggOmitFNlist),
                countsDict['badLines'],  
                countsDict['retained'], config.cleanedIsoggFN, numAltNames, 
                countsDict['unique'], self.config.uniqueIsoggFN))

    def writeUniqueSNPtable(self):
        'sort unique snp list by phylogeny and position; write to file'
        uniqueIsoggFN = self.config.uniqueIsoggFN
        
        self.snpList = sorted(self.snpList, key = attrgetter('position'))
        self.snpList = sorted(self.snpList, key = methodcaller('getDFSrank'))
        with open(uniqueIsoggFN, 'w') as uniqueIsoggFile:
            for snp in self.snpList:
                uniqueIsoggFile.write('%s\n' % snp.strWithAllNames())
    
    def isoggRecordIsBad(self, name, haplogroup, position, mutation):
        'identifies glitches with ISOGG database'
        
        if name.endswith('^'):
            self.isoggCountsDict['qc'] += 1
            return True
        
        if haplogroup.find('~') >= 0:
            self.isoggCountsDict['approxLoc'] += 1
            return True
    
        if    haplogroup.find('Investigation') >= 0 \
           or haplogroup.find('Notes') >= 0 \
           or haplogroup.find('Private') >= 0 \
           or haplogroup.find('Removed') >= 0 \
           or haplogroup.find('Withdrawn') >= 0 \
           or haplogroup.find('Freq. Mut.') >= 0 \
           or len(haplogroup) < 1:
            self.isoggCountsDict['provisional'] += 1
            return True
    
        if    len(mutation) != 4 \
           or mutation.find('?') >= 0 \
           or position.find('..') >= 0:
            self.isoggCountsDict['nonSNP'] += 1
            return True
        
        if int(position) in self.multiAllelicOldPosSet:
            self.isoggCountsDict['multiallelic'] += 1
            return True
        
        if name in self.snpNameSet:
            self.isoggCountsDict['duplicatedNames'] += 1
            return True
        
        if (position, mutation) in self.isoggOmitSet:
            self.isoggCountsDict['omitted'] += 1
            return True
    
        try:
            position = int(position)
        except ValueError:
            self.errAndLog('\nERROR. Invalid position: %s\n' % position)
            return True
        
        return False

    def checkMultiAllelics(self):
        if len(self.multiAllelicNewPosSet) > 0:
            numMulti = len(self.multiAllelicNewPosSet)
            with open(self.config.multiAllelicFoundFN, 'w') as outFile:
                for position in sorted(list(self.multiAllelicNewPosSet)):
                    outFile.write('%8d\n'% position)
            self.errAndLog('%s*** Dectected %d multiallelic positions. ***\n\n' % \
                                (utils.DASHES, numMulti) + \
                           'Please do the following and then re-run:\n' + \
                           '    cat %s >> %s\n\n' % \
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
        
        for numCharsToChop in xrange(1, len(haplogroup)):
            ancestorString = haplogroup[:-numCharsToChop]
            if ancestorString in self.hg2nodeDict:
                ancestor = self.hg2nodeDict[ancestorString]
                return ancestor.serialSplit(haplogroup)
        
        self.errAndLog('Unplaceable haplogroup: %s' % haplogroup)
