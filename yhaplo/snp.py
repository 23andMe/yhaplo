# -*- coding: utf-8 -*-
# David Poznik
# 2015.12.29
# snp.py
#
# Defines two classes:
# - SNP
# - PlatformSNP (a subclass of SNP)
#----------------------------------------------------------------------
from __future__ import absolute_import
import re
import six
import sys
from operator import attrgetter
from six.moves import range

from . import utils

class SNP(object):
    '''
    A snp knows its:
    - names
    - haplogroup
    - physical info: position, ancestral, derived
    '''

    tree = None
    config = None
    args = None
    errAndLog = None

    def __init__(self, name, haplogroup, position, ancestral, derived):
        if SNP.tree is None:
            sys.exit('ERROR. Before instantiating, must call SNP.setClassVariables(tree).')

        self.setLabel(name)
        self.nameList = [name]
        self.isRepresentative = name in SNP.tree.representativeSNPnameSet

        self.haplogroup = haplogroup
        self.position = position
        self.ancestral = ancestral
        self.derived = derived
        self.alleleSet = {ancestral, derived}

        self.node = SNP.tree.findOrCreateNode(haplogroup)
        self.node.addSNP(self)

    def setLabel(self, label):
        'sets the label and associated ivars'

        self.label = label
        self.labelLettersRank, self.labelLetters, self.labelNumber = (
            SNP.parseLabel(label, SNP.config.snpLabelLettersRankDict))
        self.labelCleaned = SNP.cleanLabel(label)

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

    @property
    def DFSrank(self):
        'returns depth-first search rank'

        return self.node.DFSrank

    @property
    def hgSNP(self):
        'string representation: truncated haplogroup label with SNP label. e.g., R-V88'

        return '%s-%s' % (self.node.hgTrunc, self.labelCleaned)

    def isDerived(self, geno):
        return geno == self.derived

    def isAncestral(self, geno):
        return geno == self.ancestral

    def isOnPlatform(self, platformVersion):
        return self.position in PlatformSNP.platformPosSetDict[platformVersion]

    def backTracePath(self):
        'returns the backtrace path (node list) for the corresponding node'

        return self.node.backTracePath()

    def addName(self, name):
        'adds name to list and updates label if appropriate'

        self.nameList.append(name)
        if name in SNP.tree.representativeSNPnameSet:
            self.isRepresentative = True

        if SNP.isApreferredName(self.label):
            if SNP.isApreferredName(name):
                SNP.errAndLog('WARNING. Two preferred names for one SNP: ' +
                              '%s, %s\n' % (name, self.label))
        elif SNP.isApreferredName(name):
            self.setLabel(name)
        else:
            labelLettersRank, labelLetters, labelNumber = (
                SNP.parseLabel(name, SNP.config.snpLabelLettersRankDict))
            if (labelLettersRank < self.labelLettersRank
                or (labelLetters == self.labelLetters
                    and labelNumber < self.labelNumber)):
                self.setLabel(name)

    @staticmethod
    def isApreferredName(name):
        '''checks wither a SNP name is in the set of preferred names,
            with or without an extension (e.g., '.1') if present'''

        return (name in SNP.tree.preferredSNPnameSet
                or name.split('.')[0] in SNP.tree.preferredSNPnameSet)

    @staticmethod
    def cleanLabel(label):
        'removes superfluous text and hyphens from a SNP label'

        for superfluousSNPtext in SNP.config.superfluousSNPtextList:
            label = label.replace(superfluousSNPtext, '')

        label = six.ensure_text(label)
        label = label.replace('-', '_').replace('^', '')
        label = six.ensure_text(label.replace(u'≤', '<=').encode('utf-8'))

        return label

    @staticmethod
    def parseLabel(name, snpLabelLettersRankDict):
        '''
        returns the priority rank of a snp name
        and a decomposition of the name into letters and a number
        '''

        match = re.search(r'([a-zA-Z-]*)([0-9]*)', str(name))
        labelLetters, labelNumber = match.group(1), match.group(2)
        labelNumber = int(labelNumber) if len(labelNumber) > 0 else 0
        if labelLetters in snpLabelLettersRankDict:
            labelLettersRank = snpLabelLettersRankDict[labelLetters]
        else:
            labelLettersRank = len(snpLabelLettersRankDict)  # max value

        return labelLettersRank, labelLetters, labelNumber

    @staticmethod
    def setClassVariables(tree):
        'enables SNP class to know about the tree instance, config, and args'

        SNP.tree = tree
        SNP.config = tree.config
        SNP.args = tree.args
        SNP.errAndLog = tree.config.errAndLog

        if SNP.config.runFromAblocks or SNP.args.writePlatformTrees:
            PlatformSNP.buildPlatformPosSetDict()
        if SNP.config.runFromAblocks:
            PlatformSNP.buildPlatformSNPlistDict()

    @staticmethod
    def prioritySortMarkerList(markerList):
        '''
        sorts a list of markers by priority ranking, with preference given to
        those deemed representative for the corresponding haplogroup
        '''

        markerList = sorted(markerList,
            key=attrgetter('labelLettersRank', 'labelLetters', 'labelNumber'))
        markerList = sorted(markerList,
            key=attrgetter('isRepresentative'), reverse=True)

        return markerList

    @staticmethod
    def mostHighlyRankedMarkerOnList(markerList):
        '''
        returns the most highly ranked marker on a list.
        the purpose of this method is to record the fact that marker lists are
        sorted with highest priority first
        '''

        if markerList:
            return markerList[0]
        else:
            return None

#--------------------------------------------------------------------------

class PlatformSNP(object):
    'A platform SNP knows its: position and ablock index'

    platformPosSetDict = dict()
    platformSNPlistDict = dict()

    def __init__(self, position):
        self.position = position
        self.ablockIndexList = SNP.config.pos2ablockIndexListDict[position]

    def __str__(self):
        return '%8d %7d' % (self.position, self.ablockIndex)

    def getConsensusGenotypeFromAblock(self, ablock):
        'given an ablock, returns a consensus genotype for this SNP'

        genotypeList = [PlatformSNP.getGenotypeFromAblock(ablock, ablockIndex)
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
        gets genotype from ablock

        input:  ablock      : a numpy array of {0, ..., 15}
                ablockIndex
        output: genotype
        '''

        if ablockIndex < len(ablock):
            diploidGenotype = SNP.config.ablockCodeToGenotypeDict[ablock[ablockIndex]]
            if diploidGenotype in SNP.config.homozygousGenotypeSet:
                return diploidGenotype[0]

        return SNP.config.missingGenotype

    @staticmethod
    def buildPlatformPosSetDict():
        'reads files to build a dictionary: platform_version -> set of positions'

        SNP.errAndLog('%sReading platform positions...\n\n' % utils.DASHES)
        for platform_version in range(1, SNP.config.maxPlatformVersionPlusOne):
            platform_pos_fn = SNP.config.platform_pos_fn_tp.format(platform_version)
            platform_pos_set = utils.readPositionsSet(
                platform_pos_fn, logFunction = SNP.errAndLog)
            PlatformSNP.platformPosSetDict[platform_version] = platform_pos_set

        SNP.errAndLog('\n')

    @staticmethod
    def buildPlatformSNPlistDict():
        'builds dictionary of platformSNP lists. key = platformVersion'

        SNP.errAndLog('Building dictionary of platform SNP lists...\n\n')
        for platformVersion in range(1, SNP.config.maxPlatformVersionPlusOne):
            platformPosSet = PlatformSNP.platformPosSetDict[platformVersion]
            platformSNPlist = list()
            for position in platformPosSet:
                platformSNPlist.append(PlatformSNP(position))

            PlatformSNP.platformSNPlistDict[platformVersion] = platformSNPlist

#--------------------------------------------------------------------------

class DroppedMarker(object):
    '''
    a marker not used for classification but potentially useful for node labeling
    examples: non-SNPs, multiallelic SNPs, and SNPs not meeting ISOGG quality guidelines
    '''

    def __init__(self, name, haplogroup):
        self.name = SNP.cleanLabel(name)
        self.haplogroup = haplogroup

    def addToNode(self):
        'adds this dropped marker to the corresponding node, if it exists'

        if self.haplogroup in SNP.tree.hg2nodeDict:
            self.setSortVariables()
            SNP.tree.hg2nodeDict[self.haplogroup].addDroppedMarker(self)
            return True
        else:
            return False

    def setSortVariables(self):
        'set variables used for priority sorting'

        self.labelLettersRank, self.labelLetters, self.labelNumber = (
            SNP.parseLabel(self.name, SNP.config.snpLabelLettersRankDict))
        self.isRepresentative = self.name in SNP.tree.representativeSNPnameSet
