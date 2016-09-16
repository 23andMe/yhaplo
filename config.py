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
    mutiCharHgShortFormString  = 'A00 A0-T A1a A1b A1b1 E1a E1b GHIJK HIJK IJK ' \
                                 'NO1 NO2 K2b Q1a Q1b R1a R1b K2c K2d K2e'
    snpLabelLettersRankString  = 'M P V L CTS AF B F Page PF Z SK'  # for prioritization
    superfluousSNPprefixString = 'IMS-'                    # stripped off snp names
    callingProgressEarlySet    = { 100, 500, 1000, 5000 }  # for progress messages
    callingProgressInterval    = 10000                     #     (23andMe only)
   
    # derived constants
    maxPlatformVersionPlusOne = maxPlatformVersion + 1
    newickSemanticTokenSet = set(newickSemanticTokenString)
    mutiCharHgShortFormSet = set(mutiCharHgShortFormString.split())
    superfluousSNPprefixList = set(superfluousSNPprefixString.split())

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
    # file names
    fileRoot               = os.path.dirname(os.path.realpath(__file__))
    inDir                  = os.path.join(fileRoot, 'input')
    outDir                 = os.path.join(fileRoot, 'output')   # project-specific files -> self.outDir

    # input | phylogenetic data
    primaryTreeFN          =  '%s/y.tree.primary.%s.nwk'   % (inDir, isoggDate)
    isoggFN                =  '%s/isogg.%s.txt'            % (inDir, isoggDate)
    isoggCorrectionsFNlist = ['%s/isogg.correct.coordinate.txt'      % inDir,
                              '%s/isogg.correct.polarize.txt'        % inDir]
    isoggOmitFNlist        = ['%s/isogg.omit.bad.txt'                % inDir,
                              '%s/isogg.omit.bad.23andMe.txt'        % inDir,
                              '%s/isogg.omit.branch.conflict.txt'    % inDir]
    isoggMultiAllelicFN    =  '%s/isogg.multiallelic.txt'            % inDir
    pagesFN                =  '%s/product/23andMe.content.pages.txt' % fileRoot

    # input | platform
    platformPosFNtp        = '%s/platform/output/v%%d.sites.txt' % fileRoot

    # input | test data
    thousandYdataFNtp      = '%s/1000Y/process/output/%%s' % fileRoot
    thousandYhgFN          = '%s/1000Y/haplogroups/4.haplogroups.called.txt' % fileRoot

    # output | phylogenetic info
    phyloOutputFNtps = {
        'withOutdirAndIsoggDate': {
            'alignedPrimaryTreeFN': '%s/y.tree.primary.aligned.%s.nwk',
            'treeFN': '%s/y.tree.%s.nwk',
            'alignedTreeFN': '%s/y.tree.aligned.%s.nwk',
            'platformTreeFNtp': '%s/y.tree.platform.%%d.%s.nwk',
            'bfTreeFN': '%s/y.tree.bf.traversal.%s.txt',
            'dfTreeFN': '%s/y.tree.df.traversal.%s.txt',
            'treeTableFN': '%s/y.tree.table.%s.txt',
            'bfPrimaryTreeFN': '%s/y.tree.primary.bf.traversal.%s.txt',
            'dfPrimaryTreeFN': '%s/y.tree.primary.df.traversal.%s.txt',
            'cleanedIsoggFN': '%s/isogg.snps.cleaned.%s.txt',
            'uniqueIsoggFN': '%s/isogg.snps.unique.%s.txt',
            'droppedIsoggFN': '%s/isogg.snps.dropped.%s.txt',
        },
        'withOutdir': {
            'multiAllelicFoundFN': '%s/multiallelic.pos',
            'pageMappingsFN': '%s/23andMe.content.page.mappings.txt',
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

    def __init__(self, description, outDir=None):
        self.args = setCommandLineArgs(description)
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

    def setDefaultAndDerivedParams(self, outDir):
        'set default and derived parameter values'
        self.outDir            = outDir if outDir is not None else Config.outDir
        self.phyloOutDir       = self.outDir
        self.vcfStartCol       = Config.vcfStartCol
        self.numCharsToCompare = Config.numCharsToCompareDefault

        # zero or one of these four will be set to True
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
        self.outDir = Config.outDir + '.1000Y'
        self.args.alternativeRoot = 'A0-T'
        self.args.importPrevCalledHaplogroup = True
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
            self.outDir = Config.outDir + '.23andMe'
        elif self.runFromVCF4:
            self.vcfStartCol = Config.vcf4startCol

    #----------------------------------------------------------------------
    def makeOutputDirectories(self):
        'make output directories if they do not already exist'
                    
        utils.mkdirP(Config.outDir)
        utils.mkdirP(self.outDir)

    #----------------------------------------------------------------------
    def setOutputFileNamesAndOpenSome(self):
        '''set log and output file names.
            open those to which we will be writing in real time'''

        for tpattr, tp in Config.phyloOutputFNtps['withOutdirAndIsoggDate'].iteritems():
            setattr(self, tpattr, tp % (self.phyloOutDir, self.isoggDate))

        for tpattr, tp in Config.phyloOutputFNtps['withOutdir'].iteritems():
            setattr(self, tpattr, tp % self.phyloOutDir)

        if self.args.singleSampleID:
            self.outFNlabel = '%s%s.' % (self.outFNlabel, self.args.singleSampleID)
        self.haplogroupCallsFN = self.constructFileName(Config.haplogroupCallsFNtp)
        self.logFN             = self.constructFileName(Config.logFNtp)

        if self.args.writeAncDerCounts:
            self.countsAncDerFN = self.constructFileName(Config.countsAncDerFNtp)
        if self.args.writeHaplogroupPaths:
            self.haplogroupPathsFN = self.constructFileName(Config.haplogroupPathsFNtp)
        if self.args.writeDerSNPs:
            self.derSNPsFN = self.constructFileName(Config.derSNPsFNtp)
        if self.args.writeDerSNPsDetail:
            self.derSNPsDetailFN = self.constructFileName(Config.derSNPsDetailFNtp)
        if self.args.writeAncSNPs:
            self.ancSNPsFN = self.constructFileName(Config.ancSNPsFNtp)
        if self.args.writeAncSNPsDetail:
            self.ancSNPsDetailFN = self.constructFileName(Config.ancSNPsDetailFNtp)
            
        # files written to in real time. open now.
        if self.args.writeHaplogroupsRealTime:
            self.haplogroupRealTimeFN = self.constructFileName(Config.haplogroupRealTimeFNtp)
            self.haplogroupRealTimeFile = open(self.haplogroupRealTimeFN, 'w', 0)
        
        if self.args.haplogroupToListGenotypesFor:
            self.hgGenosFN = self.constructFileName(Config.hgGenosFNtp) % \
                                            self.args.haplogroupToListGenotypesFor
            self.hgGenosFile = open(self.hgGenosFN, 'w', 0)

    def constructFileName(self, FNtp):
        'returns a filename, given a template, directory, and label'
        
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
        
        self.logFile = open(self.logFN, 'w', 0)
        self.errAndLog('\n%s   Y-chromosome haplogroup caller\n' % utils.DASHES + \
                       '      Command: %s\n' % ' '.join(sys.argv) + \
                       '      Log:     %s\n%s' % (self.logFN, utils.DASHES))

    def errAndLog(self, message):
        'output a message to stderr and write to the log file'
        
        sys.stderr.write(message)
        if hasattr(self, 'logFile'):
            self.logFile.write(message)
        else:
            sys.exit('ERROR. Log file not opened.')

    #----------------------------------------------------------------------
    def set23andMeArgs(self):
        'set arguments for 23andMe data'
        
        self.numCharsToCompare = 1
        self.noAblocksFN   = self.constructFileName(Config.noAblocksFNtp)
        self.noGenotypesFN = self.constructFileName(Config.noGenotypesFNtp)
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
        '''builds a dictionary with:
            key:   physical coordinate
            value: list of ablock indexes'''
        
        snpMetaDS = snpMetadata.PythonMetaDataFactory().get_dataset(Config.snpMetaDSname)
        snpMetaColList   = Config.snpMetaColList
        snpMetaArrayDict = snpMetaDS.load(snpMetaColList)
        
        self.pos2ablockIndexListDict = defaultdict(list)
        for ablockIndex, (chromosome, position) in \
                enumerate(zip(*[snpMetaArrayDict[column] for column in snpMetaColList])):
            if chromosome == Config.chromosomeInteger:
                self.pos2ablockIndexListDict[position].append(ablockIndex)
                
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
            
        self.logFile.close()

def setCommandLineArgs(description):
    parser = argparse.ArgumentParser(description=description, 
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
    group.add_argument('-pt', '--platformTrees', 
        dest='writePlatformTrees', action='store_true', default=False, 
        help='write trees whose branch lengths are '
                     'numbers of platform sites')
    group.add_argument('-m', '--mrca', type=str, nargs=2, 
        dest='mrcaHaplogroupList', 
        metavar=('haplogroup1', 'haplogroup2'),
        help='output mrca of two haplogroups\n')
    group.add_argument('-cm', '--contentMapping', 
        dest='writeContentMappings', action='store_true', default=False, 
        help='assign to each node the most closely related '
             '23andMe content page')
    
    groupDescription = 'write to files details of haplogroup calling for each sample\n' \
        '    end: written after all haplogroups have been called and samples sorted\n' \
        '    real time: written as haplogroups are called'
    group = parser.add_argument_group('details', groupDescription)
    group.add_argument('-c', '--ancDerCounts', 
        dest='writeAncDerCounts', action='store_true', default=False, 
        help='end: counts of ancestral and derived alleles encountered\n'
             'at each node visited (omits nodes with zero of each)')
    group.add_argument('-hp', '--haplogroupPaths', 
        dest='writeHaplogroupPaths', action='store_true', default=False,
        help='end: sequence of branch labels from root to call,\n'
             'with counts of derived SNPs observed')
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
             'so can this file can be sorted ex post facto with: sort -nk5')
    group.add_argument('-hg', '--hgGenos', type=str, 
        dest='haplogroupToListGenotypesFor', metavar='hg', default=None,
        help='real time: genotypes at SNPs associated with this haplogroup\n'
             'if the corresponding node is visited')
    
    groupDescription = 'restrict input or traversal'
    group = parser.add_argument_group('restrictions', groupDescription)
    group.add_argument('-po', '--primaryOnly', 
        action='store_true', default=False,
        help='do NOT import ISOGG SNPs')
    group.add_argument('-r', '--root', type=str,
        dest='alternativeRoot', metavar='hg',
        help='start searching tree from this branch')
    group.add_argument('-s', '--singleSample', type=str, 
        dest='singleSampleID', metavar='ID', 
        help='restrict to a single sample (resid for 23andMe data)\n')
    
    groupDescription = 'other options'
    group = parser.add_argument_group('other', groupDescription)
    group.add_argument('-ih', '--importPrevCalledHaplogroup', 
        action='store_true', default=False,
        help='import previously called haplogroups for comparison')
    group.add_argument('-fn', '--fileNamesOnly', 
        action='store_true', default=False,
        help='print file names to stdout and exit\n\n')
    
    return parser.parse_args()
