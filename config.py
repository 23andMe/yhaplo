# David Poznik
# 2016.01.12
# config.py
#
# Defines the Config class, which includes command-line arguments.
#----------------------------------------------------------------------
import argparse
import os
import sys
from collections import namedtuple, defaultdict

import utils

DESCRIPTION = '''
This software identifies the Y-chromosome haplogroup of each male in a sample of 
one to millions. Sequence data will yield the most highly resolved classifications, 
but the algorithm also works well with chip-based genotype data, provided a reasonable 
number of phylogenetically informative sites have been assayed. 
'''

class Config(object):
    'container for parameters, constants, filenames, and command-line arguments'

    # parameters and constants
    isoggDate       = '2016.01.04'  # date ISOGG website scraped to isoggFN
    rootHaplogroup           = 'A'  # haplogroup to associate with root node
    ancStopThresh            =   2  # BFS stopping condition parameter
    derCollapseThresh        =   2  # BFS collapsing parameter
    maxPlatformVersion       =   4  # highest 23andMe chip version
    chromosomeInteger        =  24  # for 23andMe SNP metadata
    missingGenotype          = '.'  # for text input
    missingHaplogroup        = '.'  # for output
    vcfStartCol              =   9  # first data column in .vcf
    vcf4startCol             =   7  # first data column in ".vcf4"
    numCharsToCompareDefault =   3  # for matchFlag in Sample.__str__
    
    newickSemanticTokenString  = '(),:;'                   # used in regex
    allelesString              = 'A C G T D I'             # space of possible alleles
    multiCharHgTruncString     = 'A00 A0-T A0 A1 A1a A1b A1b1 ' \
                                 'BT CT DE CF GHIJK HIJK IJK IJ LT NO'
    snpLabelLettersRankString  = 'M P V L CTS AF B F Page PF Z SK'  # for prioritization
    superfluousSNPtextList     = ['IMS-', '-null']         # stripped out of snp names
    ttamHgCallReplacementDict  = { 'BT': 'B' }             # prevents artifactual calls
    callingProgressEarlySet    = { 100, 500, 1000, 5000 }  # for progress messages
    callingProgressInterval    = 10000                     #     (23andMe only)
   
    # derived constants
    maxPlatformVersionPlusOne = maxPlatformVersion + 1
    newickSemanticTokenSet    = set(newickSemanticTokenString)
    multiCharHgTruncSet       = set(multiCharHgTruncString.split())
    multiCharHgTruncMaxLen    = max([len(elem) for elem in multiCharHgTruncSet])
    superfluousSNPtextList    = set(superfluousSNPtextList)

    alleleSet = set(allelesString.split())
    homozygousGenotypeSet = set()
    for allele in alleleSet:
        homozygousGenotypeSet.add('%s%s' % (allele, allele))

    snpLabelLettersRankDict = dict()
    for rank, letters in enumerate(snpLabelLettersRankString.split()):
        snpLabelLettersRankDict[letters] = rank

    # 23andMe dataset and feature names. Col = Column = Feature
    customerMetaDSname   = 'Metadata.customer'
    snpMetaDSname        = 'master_v4'
    ablockDSname         = 'genotype.customer.ablock'
    platformMasterDSname = 'genotype.platform:MASTER'
    
    platformVerColList   = ['is_v%d' % platformVersion \
                            for platformVersion in xrange(1, maxPlatformVersionPlusOne)]
    customerMetaColList  = ['resid', 'y_haplogroup'] + platformVerColList
    customerSexColList   = ['sex', 'sex_x', 'sex_y']
    snpMetaColList       = ['platform:chrom', 'platform:pos']
    
    # structs
    CustomerTuple = namedtuple('CustomerTuple', customerMetaColList)


    #----------------------------------------------------------------------
    # files
    
    # directories
    fileRoot      = os.path.dirname(os.path.realpath(__file__))
    inDir         = os.path.join(fileRoot, 'input')
    defaultOutDir = 'output'

    # input | phylogenetic data
    primaryTreeFN          =  '%s/y.tree.primary.%s.nwk'                  % (inDir, isoggDate)
    isoggFN                =  '%s/isogg.%s.txt'                           % (inDir, isoggDate)
    isoggCorrectionsFNlist = ['%s/isogg.correct.coordinate.txt'           % inDir,
                              '%s/isogg.correct.polarize.txt'             % inDir]
    isoggOmitFNlist        = ['%s/isogg.omit.bad.txt'                     % inDir,
                              '%s/isogg.omit.bad.23andMe.txt'             % inDir,
                              '%s/isogg.omit.branch.conflict.txt'         % inDir]
    isoggMultiAllelicFN    =  '%s/isogg.multiallelic.txt'                 % inDir
    isoggRepSNPfn          =  '%s/representative.SNPs.isogg.2015tree.txt' % inDir
    otherRepSNPfn          =  '%s/representative.SNPs.additional.txt'     % inDir
    pagesFN                =  '%s/23andMe.content.pages.txt'              % inDir

    # input | platform
    platformPosFNtp        = '%s/platform/output/v%%d.sites.txt' % fileRoot

    # input | test data
    thousandYdataFNtp      = '%s/1000Y/process/output/%%s' % fileRoot
    thousandYhgFN          = '%s/1000Y/haplogroups/4.haplogroups.called.txt' % fileRoot

    # output | phylogenetic info
    phyloOutputFNtpDictDict = {
        'withOutdirAndIsoggDate': {
            'alignedPrimaryTreeFN': '%s/y.tree.primary.aligned.%s.nwk',
            'treeFN':               '%s/y.tree.%s.nwk',
            'alignedTreeFN':        '%s/y.tree.aligned.%s.nwk',
            'platformTreeFNtp':     '%s/y.tree.platform.%%d.%s.nwk',
            'bfTreeFN':             '%s/y.tree.bf.traversal.%s.txt',
            'dfTreeFN':             '%s/y.tree.df.traversal.%s.txt',
            'treeTableFN':          '%s/y.tree.table.%s.tsv',
            'bfPrimaryTreeFN':      '%s/y.tree.primary.bf.traversal.%s.txt',
            'dfPrimaryTreeFN':      '%s/y.tree.primary.df.traversal.%s.txt',
            'cleanedIsoggFN':       '%s/isogg.snps.cleaned.%s.txt',
            'uniqueIsoggFN':        '%s/isogg.snps.unique.%s.txt',
            'droppedIsoggFN':       '%s/isogg.snps.dropped.%s.txt',
        },
        'withOutdir': {
            'multiAllelicFoundFN':  '%s/multiallelic.pos',
            'pageMappingsFN':       '%s/23andMe.content.page.mappings.txt',
            'pageUpdatesFN':        '%s/23andMe.page.updates.txt',
        }
    }

    # output | haplogroup calls, log, optional files, 23andMe auxiliary files
    logFNtp                = '%s/log.%stxt'
    haplogroupCallsFNtp    = '%s/haplogroups.%stxt'
    haplogroupRealTimeFNtp = '%s/haplogroups.realTime.%stxt'
    countsAncDerFNtp       = '%s/counts.ancDer.%stxt'
    haplogroupPathsFNtp    = '%s/paths.%stxt'
    derSNPsFNtp            = '%s/derived.snps.%stxt'
    derSNPsDetailFNtp      = '%s/derived.snps.detail.%stxt'
    ancSNPsFNtp            = '%s/ancestral.snps.%stxt'
    ancSNPsDetailFNtp      = '%s/ancestral.snps.detail.%stxt'
    hgGenosFNtp            = '%s/hg.%%s.genotypes.%stxt'
    noAblocksFNtp          = '%s/ignored.noAblocks.%sresid.txt'
    noGenotypesFNtp        = '%s/ignored.noGenotypes.%sresid.txt'

    def __init__(self, outDir=None, useDefaultCmdLineArgs=False, suppressOutputAndLog=False):
        self.useDefaultCmdLineArgs = useDefaultCmdLineArgs
        self.suppressOutputAndLog = suppressOutputAndLog
        
        self.setCommandLineArgs()
        self.setDefaultAndDerivedParams(outDir)
        self.setParamsBasedOnInputType()
        self.setParamsBasedOnRunType()
        self.makeOutputDirectories()
        self.setOutputFileNamesAndOpenSome()
        
        if self.args.fileNamesOnly:
            self.printFileNamesAndExit()
            
        self.openLogAndWelcome()
        
        if self.runFromAblocks:
            self.set23andMeArgs()
            self.get23andMeDatasets()
        
        if self.suppressOutputAndLog:
            self.overrideOutputGeneratingArgs()
    
    #----------------------------------------------------------------------
    def setDefaultAndDerivedParams(self, outDir):
        'set default and derived parameter values'
        
        # output directories
        if outDir is not None:
            self.outDir = outDir
        elif self.args.outDir is not None:
            self.outDir = self.args.outDir
        else:
            self.outDir = Config.defaultOutDir
        self.phyloOutDir = self.outDir
        
        # run type: zero or one of these four will be set to True
        self.runFromAblocks        = False
        self.runFromSampleMajorTxt = False
        self.runFromVCF            = False
        self.runFromVCF4           = False

        # derived parameters
        self.test1000Y = self.args.test1000Yall or \
                         self.args.test1000YplatformVersion or \
                         self.args.test1000Ysubset or \
                         self.args.test1000YoneID
        self.test1000YformatSpecified = self.args.test1000Yvcf or \
                                        self.args.test1000Yvcf4
                                        
        # other parameters
        self.vcfStartCol        = Config.vcfStartCol
        self.numCharsToCompare  = Config.numCharsToCompareDefault
        self.prevCalledHgFN     = self.args.prevCalledHgFN
        self.compareToPrevCalls = self.prevCalledHgFN is not None \
                                  or self.test1000Y \
                                  or self.args.compareToMetadata
        
    #----------------------------------------------------------------------
    def setParamsBasedOnInputType(self):
        'set arguments based on input type'
        
        if self.test1000Y:
            self.set1000Yargs()
        elif self.test1000YformatSpecified:
            sys.exit('ERROR. only specify test-data format when running on test data.')
        elif self.args.dataFN:
            self.parseDataFN()
        elif self.args.ablocks:
            self.runFromAblocks = True
            self.outFNlabel = '23andMe.all.'
        else:
            self.outFNlabel = ''
            
    def set1000Yargs(self):
        'set arguments for 1000Y testing'

        # args consistency checks        
        if self.args.dataFN:
            sys.stderr.write('WARNING. Ignoring explicitly specified data file.\n')
        if self.args.ablocks:
            sys.stderr.write('WARNING. Ignoring option to run on 23andMe ablocks.\n')
        if self.args.primaryOnly:
            sys.exit('ERROR. Specify --primaryOnly XOR genotype data.\n')

        # parameter values that always hold for 1000Y
        self.outDir = self.outDir + '.1000Y'
        self.args.alternativeRoot = 'A0-T'
        if not self.prevCalledHgFN:
            self.prevCalledHgFN = Config.thousandYhgFN

        # output-file labels
        if self.args.test1000Yall:
            self.outFNlabel = '1000Y.all.'
        elif self.args.test1000Ysubset:
            self.outFNlabel = '1000Y.subset.'
        elif self.args.test1000YoneID:
            self.outFNlabel = '1000Y.%s.' % self.args.test1000YoneID
            if self.args.singleSampleID and \
               self.args.singleSampleID != self.args.test1000YoneID:
                sys.exit('ERROR. Contradiction. %s vs. %s' % \
                         (self.args.singleSampleID, self.args.test1000YoneID))
        elif self.args.test1000YplatformVersion:
            self.outFNlabel = '1000Y.all.v%d.' % self.args.test1000YplatformVersion
        
        # file name for input genotypes
        if self.args.test1000Yvcf:
            self.runFromVCF = True
            dataFNlabel = self.outFNlabel + 'vcf.gz'
        elif self.args.test1000Yvcf4:
            self.runFromVCF4 = True
            dataFNlabel = self.outFNlabel + 'vcf4'
        else:
            self.runFromSampleMajorTxt = True
            dataFNlabel = self.outFNlabel + 'genos.txt'
             
        self.args.dataFN = Config.thousandYdataFNtp % dataFNlabel
        
    def parseDataFN(self):
        'set parameters based on data file name'

        dataFN = self.args.dataFN
        if dataFN.endswith('.resid.txt'):
            self.runFromAblocks = True
            self.outFNlabel = '23andMe.%s' % utils.basenameNoEnding(dataFN, '.resid.txt')
        elif dataFN.endswith('.genos.txt'):
            self.runFromSampleMajorTxt = True
            self.outFNlabel = utils.basenameNoEnding(dataFN, '.genos.txt')
        elif dataFN.endswith('.vcf'):
            self.runFromVCF = True
            self.outFNlabel = utils.basenameNoEnding(dataFN, '.vcf')
        elif dataFN.endswith('.vcf.gz'):
            self.runFromVCF = True
            self.outFNlabel = utils.basenameNoEnding(dataFN, '.vcf.gz')
        elif dataFN.endswith('.vcf4'):
            self.runFromVCF4 = True
            self.outFNlabel = utils.basenameNoEnding(dataFN, '.vcf4')
        else:
            sys.exit('\nERROR. Unknown data type: %s\n\n' % dataFN)
            
    #----------------------------------------------------------------------
    def setParamsBasedOnRunType(self):
        'set parameters based on run type'
        
        if self.runFromAblocks:
            self.outDir = self.outDir + '.23andMe'
        elif self.runFromVCF4:
            self.vcfStartCol = Config.vcf4startCol

    #----------------------------------------------------------------------
    def makeOutputDirectories(self):
        'make output directories if they do not already exist'
        
        if not self.suppressOutputAndLog:
            utils.mkdirP(self.outDir)
            utils.mkdirP(self.phyloOutDir)

    #----------------------------------------------------------------------
    def setOutputFileNamesAndOpenSome(self):
        '''
        set log and output file names.
        open those to which we will be writing in real time.
        '''

        for fn, tp in Config.phyloOutputFNtpDictDict['withOutdirAndIsoggDate'].iteritems():
            setattr(self, fn, tp % (self.phyloOutDir, self.isoggDate))

        for fn, tp in Config.phyloOutputFNtpDictDict['withOutdir'].iteritems():
            setattr(self, fn, tp % self.phyloOutDir)

        if self.args.singleSampleID:
            self.outFNlabel = '%s%s.' % (self.outFNlabel, self.args.singleSampleID)
        self.haplogroupCallsFN = self.constructOutFileName(Config.haplogroupCallsFNtp)
        self.logFN             = self.constructOutFileName(Config.logFNtp)

        if self.args.writeAncDerCounts:
            self.countsAncDerFN = self.constructOutFileName(Config.countsAncDerFNtp)
        if self.args.writeHaplogroupPaths:
            self.haplogroupPathsFN = self.constructOutFileName(Config.haplogroupPathsFNtp)
        if self.args.writeDerSNPs:
            self.derSNPsFN = self.constructOutFileName(Config.derSNPsFNtp)
        if self.args.writeDerSNPsDetail:
            self.derSNPsDetailFN = self.constructOutFileName(Config.derSNPsDetailFNtp)
        if self.args.writeAncSNPs:
            self.ancSNPsFN = self.constructOutFileName(Config.ancSNPsFNtp)
        if self.args.writeAncSNPsDetail:
            self.ancSNPsDetailFN = self.constructOutFileName(Config.ancSNPsDetailFNtp)
            
        # files written to in real time. open now.
        if self.args.writeHaplogroupsRealTime:
            self.haplogroupRealTimeFN = self.constructOutFileName(Config.haplogroupRealTimeFNtp)
            self.haplogroupRealTimeFile = open(self.haplogroupRealTimeFN, 'w', 0)
        
        if self.args.haplogroupToListGenotypesFor:
            self.hgGenosFN = self.constructOutFileName(Config.hgGenosFNtp) % \
                                            self.args.haplogroupToListGenotypesFor
            self.hgGenosFile = open(self.hgGenosFN, 'w', 0)

    def constructOutFileName(self, FNtp):
        'returns an output file name, given a template'
        
        return FNtp % (self.outDir, self.outFNlabel)
    
    #----------------------------------------------------------------------
    def printFileNamesAndExit(self):
        'prints input and output file names to stdout, then exits'
        
        print 'in:  %s' % self.args.dataFN
        print 'out: %s' % self.haplogroupCallsFN
        sys.exit()

    #----------------------------------------------------------------------
    def openLogAndWelcome(self):
        'opens log file and emits a welcome message'
        
        if self.suppressOutputAndLog:
            self.logFile = None
        else:
            self.logFile = open(self.logFN, 'w', 0)
        
        self.errAndLog('\n%s   yHaplo | Y-chromosome haplogroup caller\n' % utils.DASHES)
        if not self.useDefaultCmdLineArgs:
            self.errAndLog('      Command: %s\n' % ' '.join(sys.argv))
        if not self.suppressOutputAndLog:
            self.errAndLog('      Log:     %s\n' % self.logFN)
        self.errAndLog('%s' % utils.DASHES)

    def errAndLog(self, message):
        'output a message to stderr and write to the log file'
        
        message = message.replace(Config.fileRoot + '/', '')
        sys.stderr.write(message)
        if self.logFile:
            self.logFile.write(message)

    #----------------------------------------------------------------------
    def set23andMeArgs(self):
        'set arguments for 23andMe data'
        
        self.numCharsToCompare = 1
        self.noAblocksFN   = self.constructOutFileName(Config.noAblocksFNtp)
        self.noGenotypesFN = self.constructOutFileName(Config.noGenotypesFNtp)
        self.args.writeContentMappings = True

    def get23andMeDatasets(self):
        'get 3 datasets: snp metadata, customer metadata, genotypes (ablocks)'

        try:
            from rtk23.config import init_rtk23
            from rtk23.dataset import dataset_factory
            from rtk23.lib.coregen import VALUE_TO_CALL as ablockCodeToGenotypeDict
            from clib23 import metadata as snpMetadata
        except ImportError:
            self.exit('ERROR. Cannot import 23andMe packages: rtk23 or clib23')

        self.errAndLog('\n%sAccessing 23andMe data...\n\n' % utils.DASHES + \
                       '    Initializing rtk23 ... ')
        init_rtk23()
        self.errAndLog('Done.\n')

        self.errAndLog('    SNPs: inferring ablock indexes ... ')
        self.setPos2ablockIndexListDict(snpMetadata)
        self.errAndLog('Done.\n')

        self.errAndLog('    Customers: getting metadata...\n')
        self.customerMetaDS = dataset_factory.get_dataset(Config.customerMetaDSname)

        self.errAndLog('    Genotypes: getting ablock dataset...\n')
        self.ablockDS = dataset_factory.get_dataset(Config.ablockDSname)
        self.ablockCodeToGenotypeDict = ablockCodeToGenotypeDict

    def setPos2ablockIndexListDict(self, snpMetadata):
        '''
        builds a dictionary with:
        key:   physical coordinate
        value: list of ablock indexes
        '''
        
        snpMetaDS = snpMetadata.PythonMetaDataFactory().get_dataset(Config.snpMetaDSname)
        snpMetaColList   = Config.snpMetaColList
        snpMetaArrayDict = snpMetaDS.load(snpMetaColList)
        
        self.pos2ablockIndexListDict = defaultdict(list)
        for ablockIndex, (chromosome, position) in \
                enumerate(zip(*[snpMetaArrayDict[column] for column in snpMetaColList])):
            if chromosome == Config.chromosomeInteger:
                self.pos2ablockIndexListDict[position].append(ablockIndex)
                
    #----------------------------------------------------------------------
    def overrideOutputGeneratingArgs(self):
        'turn off all auxiliary output options'
        
        self.args.traverseBF = False
        self.args.traverseDF = False
        self.args.writeTreeTable = False
        self.args.writeContentMappings = False
        self.args.writePlatformTrees = False
    
        self.args.writeAncDerCounts = False
        self.args.writeHaplogroupPaths = False
        self.args.writeDerSNPs = False
        self.args.writeDerSNPsDetail = False
        self.args.writeAncSNPs = False
        self.args.writeAncSNPsDetail = False
        
        self.args.writeHaplogroupsRealTime = False
        self.args.haplogroupToListGenotypesFor = None

    #----------------------------------------------------------------------
    def closeFiles(self):
        'close optional real-time output files and log'
        
        if self.args.writeHaplogroupsRealTime:
            self.haplogroupRealTimeFile.close()

        if self.args.haplogroupToListGenotypesFor:
            self.hgGenosFile.close()
            self.errAndLog('Wrote genotypes at SNPs associated haplogroup %s:\n' \
                           '    %s\n\n' % (self.args.haplogroupToListGenotypesFor, 
                                           self.hgGenosFN))
        
        if self.logFile:
            self.logFile.close()

    #----------------------------------------------------------------------
    def setCommandLineArgs(self):
        'reads command-line arguments or sets defaults if self.useDefaultCmdLineArgs'
        
        parser = argparse.ArgumentParser(description=DESCRIPTION,
                                         formatter_class=argparse.RawTextHelpFormatter)
    
        # data
        group = parser.add_mutually_exclusive_group()
        group.add_argument('-i', '--input', type=str,
            dest='dataFN', metavar='fileName', default=None,
            help='specify input file in one of the following formats:\n'
                 '    .resid.txt           : list of 23andMe research IDs (col 1)\n'
                 '    .genos.txt           : sample-major text data\n'
                 '                           row 1: coordinates, col 1: sample IDs\n'
                 '    .vcf, .vcf.gz, .vcf4 : snp-major text data\n'
                 'The file format inferred from its name.')
        group.add_argument('-a', '--ablocks', 
            dest='ablocks', action='store_true', default=False,
            help='run on 23andMe ablocks\n'
                 '* the 2 data options above are mutually exclusive *\n\n')
        
        # test data
        group = parser.add_mutually_exclusive_group()
        group.add_argument('-ta', '--test1000Yall', 
            action='store_true', default=False,
            help='1000Y testing: all sites, all samples')
        group.add_argument('-ts', '--test1000Ysubset', 
            action='store_true', default=False,
            help='1000Y testing: all sites, subset of samples')
        group.add_argument('-t1', '--test1000YoneID', type=str, 
            metavar='ID',
            help='1000Y testing: all sites, one sample')
        group.add_argument('-tv', '--test1000YplatformSites', type=int, 
            dest='test1000YplatformVersion', metavar='version',
            help='1000Y testing: 23andMe sites, all samples\n'
                 '* the 4 test-data options above are mutually exclusive *\n\n')
        
        # test data format
        group = parser.add_mutually_exclusive_group()
        group.add_argument('-tvcf', '--test1000Yvcf', 
            dest='test1000Yvcf', action='store_true', default=False,
            help='1000Y testing: use .vcf.gz file rather than .genos.txt')
        group.add_argument('-tvcf4', '--test1000Yvcf4', 
            dest='test1000Yvcf4', action='store_true', default=False,
            help='1000Y testing: use .vcf4 file rather than .genos.txt\n'
                 '* the 2 test-data format options above are mutually exclusive *\n')
    
        groupDescription = 'traverse tree'
        group = parser.add_argument_group('trees', groupDescription)
        group.add_argument('-b', '--breadthFirst', 
            dest='traverseBF', action='store_true', default=False, 
            help='write bread-first traversal')
        group.add_argument('-d', '--depthFirst', 
            dest='traverseDF', action='store_true', default=False, 
            help='write depth-first (pre-order) traversal')
        group.add_argument('-dt', '--depthFirstTable', 
            dest='writeTreeTable', action='store_true', default=False, 
            help='write depth-first (pre-order) traversal table')
        group.add_argument('-m', '--mrca', type=str, nargs=2, 
            dest='mrcaHaplogroupList', 
            metavar=('haplogroup1', 'haplogroup2'),
            help='output mrca of two haplogroups')
        group.add_argument('-cm', '--contentMapping', 
            dest='writeContentMappings', action='store_true', default=False, 
            help='23andMe: map each node the most recent ancestor with an info page')
        group.add_argument('-pt', '--platformTrees', 
            dest='writePlatformTrees', action='store_true', default=False, 
            help='23andMe: write trees whose branch lengths are numbers of platform sites\n')
        
        groupDescription = 'write to files details of haplogroup calling for each sample\n' \
            '    + end       : written after all haplogroups have been called and samples sorted\n' \
            '    + real time : written as haplogroups are called'
        group = parser.add_argument_group('auxiliary output', groupDescription)
        group.add_argument('-c', '--ancDerCounts', 
            dest='writeAncDerCounts', action='store_true', default=False, 
            help='end: counts of ancestral and derived alleles encountered\n'
                 '     at each node visited (omits nodes with zero of each)')
        group.add_argument('-hp', '--haplogroupPaths', 
            dest='writeHaplogroupPaths', action='store_true', default=False,
            help='end: sequence of branch labels from root to call,\n'
                 '     with counts of derived SNPs observed')
        group.add_argument('-ds', '--derSNPs', 
            dest='writeDerSNPs', action='store_true', default=False,
            help='end: lists of derived SNPs on path')
        group.add_argument('-dsd', '--derSNPsDetail', 
            dest='writeDerSNPsDetail', action='store_true', default=False,
            help='end: detailed information about each derived SNP on path')
        group.add_argument('-as', '--ancSNPs', 
            dest='writeAncSNPs', action='store_true', default=False,
            help='end: lists of ancestral SNPs encountered in search')
        group.add_argument('-asd', '--ancSNPsDetail', 
            dest='writeAncSNPsDetail', action='store_true', default=False,
            help='end: detailed information about each ancestral SNP encountered in search\n\n')
        
        group.add_argument('-rt', '--writeRealTime',  
            dest='writeHaplogroupsRealTime', action='store_true', default=False, 
            help='real time: write haplogroups in real time. includes DFS rank,\n'
                 '           to enable ex post facto sorting: sort -nk5')
        group.add_argument('-hg', '--hgGenos', type=str, 
            dest='haplogroupToListGenotypesFor', metavar='hg', default=None,
            help='real time: genotypes at SNPs associated with this haplogroup,\n'
                 '           if the corresponding node is visited')
        
        groupDescription = 'restrict input or traversal'
        group = parser.add_argument_group('restrictions', groupDescription)
        group.add_argument('-po', '--primaryOnly', 
            action='store_true', default=False,
            help='do NOT import ISOGG SNPs')
        group.add_argument('-r', '--root', type=str,
            dest='alternativeRoot', metavar='haplogroup',
            help='start searching tree from this branch')
        group.add_argument('-s', '--singleSample', type=str, 
            dest='singleSampleID', metavar='ID', 
            help='restrict to a single sample (resid for 23andMe data)\n')
        
        groupDescription = 'previously called haplogroups'
        group.add_argument('-mdh', '--compareToMetadata',  
            dest='compareToMetadata', action='store_true', default=False, 
            help='23andMe: compare to previous calls')
        group.add_argument('-ph', '--prevCalledHgFN', type=str,
            dest='prevCalledHgFN', metavar='fileName', default=None,
            help='import previously called haplogroups:\n' + \
                 'ID in first column, Haplogroup in last\n')
    
        groupDescription = 'other options'
        group = parser.add_argument_group('other', groupDescription)
        group.add_argument('-fn', '--fileNamesOnly', 
            action='store_true', default=False,
            help='print file names to stdout and exit\n\n')
        group.add_argument('-o', '--outDir', type=str,
            dest='outDir', metavar='dirName', default=None,
            help='set output directory')
    
    
        if self.useDefaultCmdLineArgs:
            self.args = parser.parse_args([])
        else:
            self.args = parser.parse_args()  # reads from sys.argv[1:]
