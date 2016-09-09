# David Poznik
# 2016.1.26
# sample.py
# 
# Defines two classes:
# - Sample
# - Customer (a subclass of Sample)
#----------------------------------------------------------------------
import numpy as np
import sys
from collections import Counter
from itertools import izip
from operator import attrgetter

import utils
from snp import PlatformSNP

class Sample(object):
    '''A sample:
        - holds a genotype dictionary 
        - knows its haplogroup.'''
    
    tree       = None
    config     = None
    args       = None
    errAndLog  = None
    numAssigned  = 0
    numRootCalls = 0
    sampleList = list()
    prevCalledHaplogroupDict = dict()
    
    def __init__(self, ID, sampleIndex=None):
        self.ID = ID
        self.sampleIndex = sampleIndex      # for snp-major data
        self.pos2genoDict = dict()
        self.derSNPlist = None
        self.ancSNPlist = None
        self.mostDerivedSNP = None
        self.haplogroup = None
        self.hgSNP = None
        self.haplogroupDFSrank = None
        self.prevCalledHaplogroup = self.config.missingHaplogroup
        self.prevCalledHaplogroupDFSrank = 0
        self.ancDerCountTupleList = list()
        if Sample.args.importPrevCalledHaplogroup:
            self.setPrevCalledHaplogroup()
            
    # string representations and output
    #----------------------------------------------------------------------
    def __str__(self):
        'string representation gets extra information if previous calls exist'
        
        sampleString = '%-8s %-15s %-25s' % (self.ID, self.hgSNP, self.haplogroup)
        
        if Sample.args.importPrevCalledHaplogroup:
            prevHgStart = self.prevCalledHaplogroup[:Sample.config.numCharsToCompare]
            hgStart     = self.haplogroup[:Sample.config.numCharsToCompare]
            matchFlag = '.' if prevHgStart == hgStart else '*'
            sampleString = '%s %-25s %s' % (sampleString, 
                                            self.prevCalledHaplogroup, matchFlag)
        
        return sampleString
        
    def strSimple(self):
        'string representation with just ID and haplogroup'
        
        sampleString = '%-8s %-25s' % (self.ID, self.haplogroup)
        if Sample.args.importPrevCalledHaplogroup:
            sampleString = '%s %-25s' % (sampleString, self.prevCalledHaplogroup)

        return sampleString
    
    def strWithHaplogroupDFSrank(self):
        'string representation with haplogroup DFS rank tacked on'
        
        return '%s %5d' % (self, self.haplogroupDFSrank)
        
    def strPrevOnly(self):
        'string representation for before haplogroup has been called'
        
        return '%-8s %-25s' % (self.ID, self.prevCalledHaplogroup)
    
    def strSNPs(self, ancestral):
        'constructs a string representation with derived SNPs'
        
        if ancestral:
            snpList = self.ancSNPlist
        else:
            snpList = self.derSNPlist
        
        snpListString = ' '.join(snp.strShort() for snp in snpList)
        return '%-8s %-25s | %s' % (self.ID, self.haplogroup, snpListString)
    
    def strHaplogroupPath(self):
        'constructs a string representation with haplogroup path'
        
        if self.mostDerivedSNP:
            haplogroupList = list()
            for snp in self.derSNPlist:
                haplogroupList.append(snp.node.haplogroup)
            haplogroupCounter = Counter(haplogroupList)
            haplogroupPath = ' '.join([ \
                '%s:%d' % (node.haplogroup, haplogroupCounter[node.haplogroup]) \
                    for node in self.mostDerivedSNP.backTracePath() \
                    if node.haplogroup in haplogroupCounter
            ])
        else:
            haplogroupPath = ''
            
        return '%s | %s' % (self.strSimple(), haplogroupPath)

    def realTimeOutput(self):
        'generate real-time output if requested'
        
        if Sample.args.writeHaplogroupsRealTime:
            Sample.config.haplogroupRealTimeFile.write('%s\n' % self.strWithHaplogroupDFSrank())

        if Sample.args.haplogroupToListGenotypesFor:
            Sample.config.hgGenosFile.write('%s\n\n' % self.strSimple())

    # previously called haplogroups
    #----------------------------------------------------------------------
    def setPrevCalledHaplogroup(self):
        '''sets previously called haplogroup for testing/comparison.
            also sets corresponding DFS rank for sorting'''

        if not Sample.prevCalledHaplogroupDict:
            Sample.importPrevCalledHaplogroups()
        if self.ID in Sample.prevCalledHaplogroupDict:
            self.prevCalledHaplogroup = Sample.prevCalledHaplogroupDict[self.ID]
        else:
            Sample.errAndLog('WARNING. No previously called haplogroup for: %s\n' % self.ID)
            
        self.setPrevCalledHaplogroupDFSrank()
        
    @staticmethod
    def importPrevCalledHaplogroups():
        '''reads file with previously called haplogroups, 
            assuming first col = ID & last col = haplogroup'''
        
        with open(Sample.config.prevCalledHgFN, 'r') as prevCalledHgFile:
            for line in prevCalledHgFile:
                lineList = line.strip().split()
                ID, prevCalledHaplogroup = lineList[0], lineList[-1]
                Sample.prevCalledHaplogroupDict[ID] = prevCalledHaplogroup
                
    def setPrevCalledHaplogroupDFSrank(self, ignore=False):
        'sets depth-first search rank of previously called haplogroup'

        hg2nodeDict = Sample.tree.hg2nodeDict
        if not ignore and self.prevCalledHaplogroup != self.config.missingHaplogroup:
            haplogroupKey = self.prevCalledHaplogroup
            while haplogroupKey not in hg2nodeDict and len(haplogroupKey) > 0:
                haplogroupKey = haplogroupKey[:-1]
            if haplogroupKey in hg2nodeDict:
                self.prevCalledHaplogroupDFSrank = hg2nodeDict[haplogroupKey].DFSrank
        
    # mutaters
    #----------------------------------------------------------------------
    def addGeno(self, position, genotype):
        '''
        adds one value to the genotype dictionary
        if a contradiction is encountered, sets value to missing
        Note: there is no reason to call this method with a missing genotype, 
            this should be the only way missing values enter the dictionary.
        '''
        
        if position in self.pos2genoDict:
            if genotype != self.pos2genoDict[position]:
                self.pos2genoDict[position] = self.config.missingGenotype
        else:
            self.pos2genoDict[position] = genotype
            
    def appendAncDerCountTuple(self, node, numAncestral, numDerived):
        'stores results of search path'
        
        ancDerCountTuple = (node, numAncestral, numDerived)
        self.ancDerCountTupleList.append(ancDerCountTuple)

    def callHaplogroup(self):
        'finds path through tree and returns haplogroup'
        
        Sample.numAssigned   += 1
        path, self.ancSNPlist = Sample.tree.identifyPhylogeneticPath(self)
        self.derSNPlist       = path.derSNPlist
        self.mostDerivedSNP   = path.mostDerivedSNP

        if self.mostDerivedSNP:
            self.haplogroup   = self.mostDerivedSNP.haplogroup
            self.hgSNP        = self.mostDerivedSNP.getHgSNP()
        else:
            Sample.numRootCalls += 1
            self.haplogroup   = Sample.tree.root.haplogroup
            self.hgSNP        = self.haplogroup
            
        self.haplogroupDFSrank = Sample.tree.getDFSrank(self.haplogroup)
        self.realTimeOutput()
        self.freeUpMemory()

    def freeUpMemory(self):
        'free up some memory if possible'
        
        self.pos2genoDict = None
        if not (Sample.args.writeDerSNPs or \
                Sample.args.writeDerSNPsDetail or \
                Sample.args.writeHaplogroupPaths):
            self.derSNPlist = None
        if not (Sample.args.writeAncSNPs or \
                Sample.args.writeAncSNPsDetail):
            self.ancSNPlist = None
    

    #----------------------------------------------------------------------
    # Run: main entry point
    #----------------------------------------------------------------------
    @staticmethod
    def callHaplogroups(config, tree):
        'this method is the entry point and is to be called from outside.'
        
        Sample.setTreeConfigAndArgs(config, tree)
        Sample.testNumberOfRunModes(config)

        if config.runFromSampleMajorTxt:
            Sample.runFromSampleMajorTxt()
        elif config.runFromVCF or config.runFromVCF4:
            Sample.runFromVCF()
        elif config.runFromAblocks:
            Customer.runFromAblocks()
        else:
            Sample.errAndLog('%sNo input genotypes specified. Exiting.\n\n' % utils.DASHES)
            
        config.closeFiles()

    @staticmethod
    def setTreeConfigAndArgs(config, tree):
        'enables Sample class to know about the tree instance, config, and args'
        
        Sample.config    = config
        Sample.args      = config.args
        Sample.errAndLog = config.errAndLog
        Sample.tree      = tree
        
        if Sample.args.writeHaplogroupsRealTime:
            Sample.realTimeHaplogroupWritingMessage()

    @staticmethod
    def realTimeHaplogroupWritingMessage():
        'emit a message for real-time haplogroup writing'
        
        Sample.errAndLog( \
            '%sWill write haplogroups as they are called:\n' % utils.DASHES + \
            '    %s\n\n' % Sample.config.haplogroupRealTimeFN + \
            'Note: This file includes DFS rank, so it can be sorted ex post facto with:\n' + \
            '    sort -nk5 %s\n\n' % Sample.config.haplogroupRealTimeFN)

    @staticmethod
    def testNumberOfRunModes(config):
        'consistency check the number of run modes'
        
        numberOfRunModesSelected = config.runFromSampleMajorTxt + \
                                   config.runFromVCF + \
                                   config.runFromVCF4 + \
                                   config.runFromAblocks
        if numberOfRunModesSelected > 1:
            sys.exit('ERROR. Expecting no more than one run mode\n' + \
                     '    %d selected\n' % numberOfRunModesSelected)


    #----------------------------------------------------------------------
    # Run option 1: sample-major text data
    #----------------------------------------------------------------------
    @staticmethod
    def runFromSampleMajorTxt():
        'run pipeline on sample-major data'
        
        Sample.processSampleMajorTxtandCallHaplogroups()
        Sample.sortSampleList()
        Sample.writeSampleList()
        
    @staticmethod
    def processSampleMajorTxtandCallHaplogroups():
        '''
        reads in sample major data, calling haplogroup for each line.
        returns list of sample objects with genotype data purged.

        assumed format:
            row 1    = physical coordinates
            column 1 = sample ID
        '''
        
        genoFN               = Sample.args.dataFN
        genoFile, genoReader = utils.getCSVreader(genoFN, delimiter='\t')
        Sample.errAndLog('%sReading genotype data...\n    %s\n\n' % \
                         (utils.DASHES, genoFN))

        # determine relevant physical coordinates and corresponding columns        
        allPositionsList = [int(position) for position in genoReader.next()[1:]]
        columnPositionTupleList = list()
        for column, position in enumerate(allPositionsList):
            if position in Sample.tree.snpPosSet:
                columnPositionTupleList.append((column, position))
        
        # read genotypes, call haplogroups
        for genoList in genoReader:
            ID, genoList = genoList[0], genoList[1:]
            if Sample.args.singleSampleID and ID != Sample.args.singleSampleID:
                continue

            sample = Sample(ID)
            for column, position in columnPositionTupleList:
                genotype = genoList[column]
                if genotype != Sample.config.missingGenotype:
                    sample.addGeno(position, genotype)
                    
            sample.callHaplogroup()
            Sample.sampleList.append(sample)
            
        genoFile.close()


    #----------------------------------------------------------------------
    # Run option 2: snp-major txt data (VCF)
    #----------------------------------------------------------------------
    @staticmethod
    def runFromVCF():
        'run pipeline on snp-major data'
        
        Sample.loadDataFromVCF()
        if Sample.args.writeHaplogroupsRealTime or \
           Sample.args.haplogroupToListGenotypesFor:
            Sample.sortSampleList(sortByPrevHg=True)
        for sample in Sample.sampleList:
            sample.callHaplogroup()
        Sample.sortSampleList()
        Sample.writeSampleList()

    @staticmethod
    def loadDataFromVCF():
        'constructs list of sample objects, each with a genotype dictionary'
        
        vcfFN              = Sample.args.dataFN
        vcfFile, vcfReader = utils.getCSVreader(vcfFN, delimiter='\t')
        Sample.setSampleListFromVCFheader(vcfReader)

        Sample.errAndLog('%sReading genotype data...\n    %s\n\n' % \
                         (utils.DASHES, vcfFN))

        for lineList in vcfReader:
            position = int(lineList[1])
            if position in Sample.tree.snpPosSet:
                genoList = lineList[Sample.config.vcfStartCol:]
                for sample in Sample.sampleList:
                    genotype = genoList[sample.sampleIndex].split(':')[0]
                    
                    if genotype == Sample.config.missingGenotype:
                        continue
                    elif Sample.config.runFromVCF:      # as opposed to .vcf4
                        ref, alt = lineList[3:5]
                        if genotype == '0':
                            genotype = ref
                        elif genotype == '1':
                            genotype = alt
                            
                    sample.addGeno(position, genotype)

        vcfFile.close()
    
    @staticmethod
    def setSampleListFromVCFheader(vcfReader):
        '''reads a VCF header until it finds the sample IDs,
            then uses them to construct a list of sample objects'''
        
        for lineList in vcfReader:
            if lineList[0][:2] != '##':     # ignore metadata
                break
        Sample.validateVCFheader(lineList)

        idList = lineList[Sample.config.vcfStartCol:]
        for sampleIndex, ID in enumerate(idList):
            if not Sample.args.singleSampleID or ID == Sample.args.singleSampleID:
                sample = Sample(ID, sampleIndex)
                Sample.sampleList.append(sample)
    
    @staticmethod
    def validateVCFheader(lineList):
        'checks that the second element of a list is POS'
        
        col2label = lineList[1]
        if col2label != 'POS':
            sys.exit('ERROR. Invalid VCF. Expected column 2 header to be POS.\n' + \
                     '       Instead found: %s' % col2label)


    # sample list manipulation and writing
    #----------------------------------------------------------------------
    @staticmethod
    def sortSampleList(sortByPrevHg=False):
        'sorts sample list by haplogroup (previously called or current)'
        
        if sortByPrevHg and Sample.args.importPrevCalledHaplogroup:
            primarySortKey = 'prevCalledHaplogroupDFSrank'
        else:
            primarySortKey = 'haplogroupDFSrank'
        
        Sample.sampleList = sorted(Sample.sampleList, 
                                   key = attrgetter(primarySortKey, 'ID'))

    @staticmethod
    def writeSampleList():
        'writes haplogroup and other optional data for each sample'

        Sample.reportCounts()
        Sample.writeHaplogroups()

        if Sample.args.writeAncDerCounts:
            Sample.writeAncDerCounts()
            
        if Sample.args.writeHaplogroupPaths:
            Sample.writeHaplogroupPaths()
        
        if Sample.args.writeDerSNPs:
            Sample.writeSNPs()
        
        if Sample.args.writeDerSNPsDetail:
            Sample.writeSNPsDetail()
            
        if Sample.args.writeAncSNPs:
            Sample.writeSNPs(ancestral=True)

        if Sample.args.writeAncSNPsDetail:
            Sample.writeSNPsDetail(ancestral=True)
            
    @staticmethod
    def reportCounts():
        'report number assigned and number assigned to root'
        
        Sample.errAndLog( \
            '    %8d haplogroups assigned\n' % Sample.numAssigned + \
            '    %8d assigned to Root haplogroup: %s\n\n' % \
                    (Sample.numRootCalls, Sample.tree.root.haplogroup))
        
        if Sample.numRootCalls > 0:
            Sample.warnVariantsOnlyData()
            
    @staticmethod
    def warnVariantsOnlyData():
        'warning for datasets that exclude sites with no variation in the sample'
        
        Sample.errAndLog( \
            'WARNING. If the dataset does not include fixed reference sites,\n' + \
            '         re-run with alternative root (e.g., with: -r A0-T).\n\n\n')

    @staticmethod
    def writeHaplogroups():
        'writes haplogroup of each sample'
    
        with open(Sample.config.haplogroupCallsFN, 'w') as haplogroupCallsFile: 
            for sample in Sample.sampleList:
                haplogroupCallsFile.write('%s\n' % sample)
                
        Sample.errAndLog( \
            'Wrote called haplogroups:\n' + \
            '    %s\n\n' % Sample.config.haplogroupCallsFN)

    @staticmethod
    def writeAncDerCounts():
        '''writes counts of ancestral and derived alleles encountered
             at each node visited (excluding nodes with zero of each)'''
    
        with open(Sample.config.countsAncDerFN, 'w') as countsAncDerFile: 
            for sample in Sample.sampleList:
                for node, numAncestral, numDerived in sample.ancDerCountTupleList:
                    if numAncestral > 0 or numDerived > 0:
                        countsAncDerFile.write('%-8s %-20s %3d %3d\n' % 
                                               (sample.ID, node.label, 
                                                numAncestral, numDerived))

                countsAncDerFile.write('%s\n\n' % sample.strSimple())
                                                  
        Sample.errAndLog('Wrote counts of ancestral and derived alleles encountered\n' + \
                         'at each node visited (excluding nodes with zero of each):\n' + \
                         '    %s\n\n' % Sample.config.countsAncDerFN)

    @staticmethod        
    def writeHaplogroupPaths():
        'writes haplogroup path for each sample'

        with open(Sample.config.haplogroupPathsFN, 'w') as haplogroupPathsFile: 
            for sample in Sample.sampleList:
                haplogroupPathsFile.write('%s\n' % sample.strHaplogroupPath())

        Sample.errAndLog('Wrote sequences of haplogroups from root to calls,\n' + \
                         'with counts of derived SNPs observed:\n' + \
                         '    %s\n\n' % Sample.config.haplogroupPathsFN)
    
    @staticmethod
    def writeSNPs(ancestral=False):
        '''for each sample, writes list of derived SNPs on path
            or list of ancestral SNPs encountered in search'''
        
        if ancestral:
            outFN = Sample.config.ancSNPsFN
            typeOfSNPs = 'ancestral SNPs encountered in search'
        else:
            outFN = Sample.config.derSNPsFN
            typeOfSNPs = 'derived SNPs on path'
    
        with open(outFN, 'w') as outFile: 
            for sample in Sample.sampleList:
                outFile.write('%s\n' % sample.strSNPs(ancestral))

        Sample.errAndLog('Wrote lists of %s:\n    %s\n\n' % (typeOfSNPs, outFN))
        
    @staticmethod
    def writeSNPsDetail(ancestral=False):
        '''for each sample, writes detailed information about each 
            derived SNP on path or about each ancestral SNP encountered in search'''
        
        if ancestral:
            outFN = Sample.config.ancSNPsDetailFN
            typeOfSNPs = 'ancestral SNP encountered in search'
        else:
            outFN = Sample.config.derSNPsDetailFN
            typeOfSNPs = 'derived SNP on path'

        with open(outFN, 'w') as outFile: 
            for sample in Sample.sampleList:
                outFile.write('%s\n' % sample)
                snpList = sample.ancSNPlist if ancestral else sample.derSNPlist
                for snp in snpList:
                    outFile.write('%-8s %s\n' % (sample.ID, snp))
                outFile.write('\n')
                                                  
        Sample.errAndLog('Wrote detailed information about each %s:\n' % typeOfSNPs + \
                         '    %s\n\n' % outFN)


#--------------------------------------------------------------------------

class Customer(Sample):
    'a customer is a sample with 23andMe metadata and genotypes in ablocks'
    
    numNoAblock, numNoGenotypes = 0, 0
    noAblocksFile, noGenotypesFile = None, None
    
    def __init__(self, customerTuple):
        self.customerTuple = customerTuple
        Sample.__init__(self, customerTuple.resid)
        
    def setPrevCalledHaplogroup(self):
        '''called from Customer.__init__() via Sample.__init__()
            for testing/comparison, sets previously called haplogroup from 
            original 23andMe algorithm. does not set corresponding DFS rank
            since the nomenclature has changed substantially'''

        self.prevCalledHaplogroup = self.customerTuple.y_haplogroup
        
    def loadAblockAndCallHaplogroup(self, purgeGenotypes=True):
        'loads ablock, reads relevant genotypes, calls haplogroup'

        try:
            ablock = Sample.config.ablockDS.load_block(self.ID)
        except:
            Customer.noAblocksFile.write('%d\n' % self.ID)
            Customer.numNoAblock += 1
            return False
        
        if self.readAblockGenotypes(ablock):
            self.callHaplogroup()
            return True
        else:
            Customer.noGenotypesFile.write('%d\n' % self.ID)
            Customer.numNoGenotypes += 1
            return False
        
    def readAblockGenotypes(self, ablock):
        'pulls phylogenetically informative genotypes from ablock'
        
        hasGenotypes = False
        for platformVersion in xrange(1, Sample.config.maxPlatformVersionPlusOne):
            if getattr(self.customerTuple, 'is_v%d' % platformVersion):
                hasGenotypes = True
                platformSNPlist = PlatformSNP.platformSNPlistDict[platformVersion]
                for platformSNP in platformSNPlist:
                    genotype = platformSNP.getConsensusGenotypeFromAblock(ablock)
                    self.addGeno(platformSNP.position, genotype)

        return hasGenotypes


    #----------------------------------------------------------------------
    # Run option 3: sample-major 23andMe ablock data
    #----------------------------------------------------------------------
    @staticmethod
    def runFromAblocks():
        'run pipeline on sample-major 23andMe ablock data'
        
        Customer.openAuxiliaryOutputFiles()

        Sample.errAndLog('%sProcessing 23andMe customer data...\n\n' % utils.DASHES)
        residList = Customer.generateResidList()
        customerTupleList = Customer.buildCustomerTupleList(residList)

        Sample.errAndLog('\n%sCalling haplogroups...\n\n' % (utils.DASHES) + \
                         ' Progress...\n')
        for customerTuple in customerTupleList:
            customer = Customer(customerTuple)
            if customer.loadAblockAndCallHaplogroup(purgeGenotypes=True):
                Sample.sampleList.append(customer)
                Customer.emitProgress()
        
        Customer.closeAuxiliaryFilesAndReportCounts()
        Sample.sortSampleList()
        Sample.writeSampleList()

    @staticmethod
    def openAuxiliaryOutputFiles():
        'open auxiliary output files'
        
        Customer.noAblocksFile   = open(Sample.config.noAblocksFN,   'w')
        Customer.noGenotypesFile = open(Sample.config.noGenotypesFN, 'w')

    @staticmethod
    def generateResidList():
        '''
        3 possibilities:
            -a -s RESID           -> a single research ID has been specified
            -i FILENAME.resid.txt -> read research IDs from file
            -a                    -> return empty list to indicate no subsetting
        '''
        
        residList = list()

        if Sample.args.singleSampleID:
            resid = Customer.generateResid(Sample.args.singleSampleID)
            residList = [resid]
            Sample.errAndLog('Will call haplogroup for:\n    %d\n\n' % resid)
        elif Sample.args.dataFN:
            residFN = Sample.args.dataFN
            utils.checkFileExistence(residFN)
            Sample.errAndLog('Reading research IDs:\n    %s\n' % residFN)
            with open(residFN, 'r') as residFile:
                for lineList in residFile:
                    ID = lineList.strip().split()[0]
                    residList.append(Customer.generateResid(ID))
            Sample.errAndLog('    %8d read\n'     % len(residList))
            Sample.errAndLog('    %8d unique\n\n' % len(set(residList)))

        return residList
    
    @staticmethod
    def generateResid(ID):
        'converts ID to integer, exiting gracefully if not possible'
        
        try: 
            resid = int(ID)
        except ValueError: 
            sys.exit('\nERROR. Cannot convert ID to integer: %s' % ID)
        
        return resid
        
    @staticmethod
    def buildCustomerTupleList(residList):
        '''1. builds a customer mask.
              a. if resids have been specified, use those.
              b. otherwise, use all males.
           2. generates a list of customer tuples subsetted to this mask'''
        
        Sample.errAndLog('Building customer mask ... ')

        customerMetaDS        = Sample.config.customerMetaDS
        customerMetaArrayDict = customerMetaDS.load(Sample.config.customerMetaColList)
        
        if len(residList) > 0:
            mask = np.in1d(customerMetaArrayDict['resid'], residList)
            maskType = 'specified'
        else:
            customerSexArrayDict = customerMetaDS.load(Sample.config.customerSexColList)
            mask = (customerSexArrayDict['sex'] == 'M') & \
                   (customerSexArrayDict['sex_x'] == 'M') & \
                   (customerSexArrayDict['sex_y'] == 'M')
            maskType = 'male'

        customerTupleList = [Sample.config.CustomerTuple(*row) \
            for row in izip(*[customerMetaArrayDict[column][mask] \
                             for column in Sample.config.customerMetaColList])]
        
        Sample.errAndLog('Done.\n' +
            '    %8d %s customers to be processed\n' % \
                (len(customerTupleList), maskType))
        return customerTupleList
    
    @staticmethod
    def emitProgress():
        'emit a message indicating how many haplogroups have been assigned thus far'
        
        if Sample.numAssigned in Sample.config.callingProgressEarlySet \
                or Sample.numAssigned % Sample.config.callingProgressInterval == 0:
            Sample.errAndLog('    %8d haplogroups assigned\n' % Sample.numAssigned)

    @staticmethod
    def closeAuxiliaryFilesAndReportCounts():
        'close auxiliary files and report counts'

        for File in [Customer.noAblocksFile, Customer.noGenotypesFile]:
            File.close()

        Sample.errAndLog('\n' + \
            '    %8d resid(s) ignored | could not load ablock:\n' % Customer.numNoAblock + \
            '             %s\n'   % Sample.config.noAblocksFN + \
            '    %8d resid(s) ignored | no genotypes:\n' % Customer.numNoGenotypes + \
            '             %s\n\n' % Sample.config.noGenotypesFN)
