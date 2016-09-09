# David Poznik
# 2015.12.29
# snp.py
# 
# Defines two classes:
# - SNP
# - PlatformSNP (a subclass of SNP)
#----------------------------------------------------------------------
import re
import sys

import utils

class SNP(object):
    '''
    A snp knows its:
        - names
        - haplogroup
        - physical info: position, ancestral, derived
    '''
    
    tree      = None
    config    = None
    args      = None
    errAndLog = None

    def __init__(self, name, haplogroup, position, ancestral, derived):
        if SNP.tree is None:
            sys.exit('ERROR. Before instantiating, must call SNP.setClassVariables(tree).')

        self.label = name
        self.labelLettersRank, self.labelLetters, self.labelNumber = \
            SNP.parseLabel(self.label, SNP.config.snpLabelLettersRankDict)
        self.nameList = [name]
        
        self.haplogroup = haplogroup
        self.position   = position
        self.ancestral  = ancestral
        self.derived    = derived
        self.alleleSet  = { ancestral, derived }
        
        self.node = SNP.tree.findOrCreateNode(haplogroup)
        self.node.addSNP(self)
        
    def __str__(self):
        'medium-length string representation'
        
        return '%-15s %-25s %8d %s->%s' % (self.label, self.node.label,
            self.position, self.ancestral, self.derived)
    
    def strWithAllNames(self):
        'long string representation: normal str plus comma-separated list of names'
        
        return '%s     %s' % (str(self), ','.join(self.nameList))
        
    def strShort(self):
        'short string representation: node label and snp label' 
        
        return '%s:%s' % (self.node.label, self.label)
        
    def isDerived(self, geno):
        return geno == self.derived
    
    def isAncestral(self, geno):
        return geno == self.ancestral
    
    def isOnPlatform(self, platformVersion):
        return self.position in PlatformSNP.platformPosSetDict[platformVersion]
        
    def addName(self, name):
        'add name to list, possibly setting label'
        
        self.nameList.append(name)
        labelLettersRank, labelLetters, labelNumber = \
            SNP.parseLabel(name, SNP.config.snpLabelLettersRankDict)
        if (labelLettersRank < self.labelLettersRank) or \
                (labelLetters == self.labelLetters and \
                 labelNumber < self.labelNumber):
            self.label = name
            self.labelLettersRank = labelLettersRank
            self.labelLetters = labelLetters
            self.labelNumber = labelNumber
            
    def getDFSrank(self):
        'returns depth-first search rank'
        
        return self.node.DFSrank
    
    def getHgSNP(self):
        'string representation with shortened haplogroup'
        
        return '%s-%s' % (self.node.hgShort, self.getTrimmedLabel())
    
    def getTrimmedLabel(self):
        'returns label with superfluous prefixes trimmed'
        
        trimmedLabel = self.label
        for superfluousSNPprefix in SNP.config.superfluousSNPprefixList:
            if trimmedLabel[:len(superfluousSNPprefix)] == superfluousSNPprefix:
                trimmedLabel = trimmedLabel[len(superfluousSNPprefix):]
            
        return trimmedLabel
    
    def backTracePath(self):
        'returns the backtrace path (node list) for the corresponding node'
        
        return self.node.backTracePath()

    @staticmethod
    def parseLabel(name, snpLabelLettersRankDict):
        '''returns the priority rank of a snp name 
            and a decomposition of the name into letters and a number'''
        
        match = re.search(r'([a-zA-Z]*)([0-9]*)', name)
        letters, number = match.group(1), match.group(2)
        number = int(number) if len(number) > 0 else 0
        if letters in snpLabelLettersRankDict:
            labelLettersRank = snpLabelLettersRankDict[letters]
        else:
            labelLettersRank = len(snpLabelLettersRankDict)  # max value
        
        return labelLettersRank, letters, number
        
    @staticmethod
    def setClassVariables(tree):
        'enables SNP class to know about the tree instance, config, and args'
        
        SNP.tree      = tree
        SNP.config    = tree.config
        SNP.args      = tree.args
        SNP.errAndLog = tree.config.errAndLog
        
        if SNP.config.runFromAblocks or SNP.args.writePlatformTrees:
            PlatformSNP.buildPlatformPosSetDict()
        if SNP.config.runFromAblocks:
            PlatformSNP.buildPlatformSNPlistDict()

#--------------------------------------------------------------------------

class PlatformSNP(object):
    'A platform SNP knows its: position and ablock index'

    platformPosSetDict  = dict()
    platformSNPlistDict = dict()
    
    def __init__(self, position):
        self.position = position
        self.ablockIndexList = SNP.config.pos2ablockIndexListDict[position]

    def __str__(self):
        return '%8d %7d' % (self.position, self.ablockIndex)
    
    def getConsensusGenotypeFromAblock(self, ablock):
        'given an ablock, returns a consensus genotype for this SNP'
        
        genotypeList = [PlatformSNP.getGenotypeFromAblock(ablock, ablockIndex) \
                        for ablockIndex in self.ablockIndexList]
        if len(genotypeList) == 1:
            genotype = genotypeList[0]
        else:
            genotypeSet = set(genotypeList)
            genotypeSet.discard(SNP.config.missingGenotype)
            if len(genotypeSet) == 1:
                genotype = genotypeSet.pop()
            else:
                genotype = SNP.config.missingGenotype
            
        return genotype

    @staticmethod
    def getGenotypeFromAblock(ablock, ablockIndex):
        '''
        input:  ablock      = a numpy array of { 0, ..., 15 }
                ablockIndex
        output: genotype
        '''
        
        diploidGenotype = SNP.config.ablockCodeToGenotypeDict[ablock[ablockIndex]]
        if diploidGenotype in SNP.config.homozygousGenotypeSet:
            return diploidGenotype[0]
        else:
            return SNP.config.missingGenotype

    @staticmethod
    def buildPlatformPosSetDict():
        'reads files to build a dictionary: platformVersion -> set of positions'
        
        SNP.errAndLog('%sReading platform positions...\n\n' % utils.DASHES)
        for platformVersion in xrange(1, SNP.config.maxPlatformVersionPlusOne):
            platformPosFN = SNP.config.platformPosFNtp % platformVersion
            platformPosSet = utils.readPositionsSet(platformPosFN, 
                                                    logFunction = SNP.errAndLog)
            PlatformSNP.platformPosSetDict[platformVersion] = platformPosSet
                        
        SNP.errAndLog('\n')

    @staticmethod
    def buildPlatformSNPlistDict():
        'builds dictionary of platformSNP lists. key = platformVersion'
        
        SNP.errAndLog('Building dictionary of platform SNP lists...\n\n')
        for platformVersion in xrange(1, SNP.config.maxPlatformVersionPlusOne):
            platformPosSet  = PlatformSNP.platformPosSetDict[platformVersion]
            platformSNPlist = list()
            for position in platformPosSet:
                platformSNPlist.append(PlatformSNP(position))
                            
            PlatformSNP.platformSNPlistDict[platformVersion] = platformSNPlist
