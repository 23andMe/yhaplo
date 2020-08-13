# David Poznik
# 2016.01.12
# config.py
#
# Defines the Config class, which includes command-line arguments.
#----------------------------------------------------------------------
from __future__ import absolute_import, print_function
import argparse
import os
import six
import sys
from collections import namedtuple, defaultdict
from six.moves import range
from six.moves import zip

from . import utils, __version__

DESCRIPTION = '''
yhaplo identifies the Y-chromosome haplogroup of each male in a sample of one
to millions. Sequence data will yield the most highly resolved classifications,
but the algorithm also works well with chip-based genotype data, provided a
reasonable number of phylogenetically informative sites have been assayed.
'''

#----------------------------------------------------------------------
# constants

ANC_STOP_THRESH_DEFAULT = 2         # BFS stopping condition parameter default
DER_COLLAPSE_THRESH_DEFAULT = 2     # BFS collapsing parameter default

#----------------------------------------------------------------------
class Config(object):
    'container for parameters, constants, filenames, and command-line arguments'

    # parameters and constants
    #--------------------------------------------------------------------
    isoggDate = '2016.01.04'        # date ISOGG website scraped to isoggFN
    rootHaplogroup = 'A'            # haplogroup to associate with root node
    missingGenotype = '.'           # for text input
    missingHaplogroup = '.'         # for output
    vcf_chrom_label_set = {'Y', 'chrY', '24'}  # to restrict VCF input to chrY
    vcfStartCol = 9                 # first data column in .vcf
    vcf4startCol = 7                # first data column in ".vcf4"
    numCharsToCompareDefault = 3    # for matchFlag in Sample.__str__

    newickSemanticTokenString = '(),:;'                             # used in regex
    allelesString = 'A C G T D I'                                   # space of possible alleles
    snpLabelLettersRankString = 'M P V L CTS AF B F Page U PF Z SK' # for prioritization
    superfluousSNPtextList = ['IMS-', '-null']                      # stripped out of snp names
    multiCharHgTruncString = 'A00 A0-T A0 A1 A1a A1b A1b1 BT CT DE CF GHIJK HIJK IJK IJ LT NO'

    # derived constants
    newickSemanticTokenSet = set(newickSemanticTokenString)
    multiCharHgTruncSet = set(multiCharHgTruncString.split())
    multiCharHgTruncMaxLen = max([len(elem) for elem in multiCharHgTruncSet])
    superfluousSNPtextList = set(superfluousSNPtextList)
    alleleSet = set(allelesString.split())
    homozygousGenotypeSet = {'%s%s' % (allele, allele) for allele in alleleSet}
    snpLabelLettersRankDict = {letters: rank
                               for rank, letters in enumerate(snpLabelLettersRankString.split())}

    # 23andMe-specific parameters and constants
    #--------------------------------------------------------------------
    ttamHgCallReplacementDict = {'BT': 'B'}                 # prevents artifactual calls
    callingProgressEarlySet = {100, 500, 1000, 5000}        # for progress messages
    callingProgressInterval = 10000                         # for progress messages

    # SNPs
    maxPlatformVersion = 5                                  # most recent chip version
    chromosomeInteger = 24
    snpMetaDSname = 'Metadata.master_v%d' % maxPlatformVersion
    snpMetaColList = ['platform:chrom', 'platform:pos']     # features to draw from SNP metadata
    maxPlatformVersionPlusOne = maxPlatformVersion + 1

    # genotypes
    ablockDSnameDefault = 'genotype.customer.ablock'
    ablock_fn_tp = '{}.ablock.npy.gz'
    bitpacked_ablock_fn_tp = '*.{}.ablock' # for use in a glob.glob

    # samples
    customerMetaDSname = 'Metadata.customer'                    # for when using ablockDSnameDefault
    customerIDcol             = 'resid'
    customerPlatformColList   = ['is_v%d' % platformVersion
                                 for platformVersion in range(1, maxPlatformVersionPlusOne)]
    customerSexColList        = ['sex', 'sex_x', 'sex_y']       # for --allMaleCustomers option
    customerPrevHaplogroupCol = 'y_haplogroup'                  # for --compareToMetadata option
    customerMetaColList = [customerIDcol] + customerPlatformColList   # for whenever metadata used
    CustomerTuple = namedtuple('CustomerTuple', customerMetaColList + [customerPrevHaplogroupCol])

    # filenames
    #----------------------------------------------------------------------
    # directories
    softwareDir = os.path.dirname(os.path.realpath(__file__))
    inDir = os.path.join(softwareDir, 'input')
    defaultOutDir = 'output'

    # input: phylogenetic data
    primaryTreeFN          =  '%s/y.tree.primary.%s.nwk'                     % (inDir, isoggDate)
    isoggFN                =  '%s/isogg.%s.txt'                              % (inDir, isoggDate)
    isoggCorrectionsFNlist = ['%s/isogg.correct.coordinate.txt'              % inDir,
                              '%s/isogg.correct.polarize.txt'                % inDir]
    isoggOmitFNlist        = ['%s/isogg.omit.bad.txt'                        % inDir,
                              '%s/isogg.omit.bad.23andMe.txt'                % inDir,
                              '%s/isogg.omit.branch.conflict.txt'            % inDir,
                              '%s/isogg.omit.branch.conflict.23andMe.v5.txt' % inDir]
    isoggMultiAllelicFN    =  '%s/isogg.multiallelic.txt'                    % inDir
    isoggRepSNPfn          =  '%s/representative.SNPs.isogg.2015tree.txt'    % inDir
    otherRepSNPfn          =  '%s/representative.SNPs.additional.txt'        % inDir
    preferredSNPnamesFN    =  '%s/preferred.snpNames.txt'                    % inDir
    pagesFN                =  '%s/23andMe.content.pages.txt'                 % inDir

    # input: platform
    platform_pos_fn_tp = os.path.join(softwareDir, 'platform_sites/v{}.sites.txt')
    
    # input: example data
    kgp_subset_fn = os.path.join(softwareDir, 'data', '1000Y.subset.genos.txt')
    kgp_single_sample_vcf_fn = os.path.join(softwareDir, 'data', 'HG01938.vcf.gz')
    
    # input: test data
    thousandYdataFNtp = os.path.join(softwareDir, '1000Y/process/output/%s')
    thousandYhgFN = os.path.join(softwareDir, '1000Y/haplogroups/4.haplogroups.called.txt')

    # output: phylogenetic info
    phyloOutputFNtpDictDict = {
        'withOutdirAndIsoggDate': {
            'alignedPrimaryTreeFN':  '%s/y.tree.primary.aligned.ycc.%s.nwk',
            'yccTreeFN':             '%s/y.tree.ycc.%s.nwk',
            'hgsnpTreeFN':           '%s/y.tree.hgSNP.%s.nwk',
            'alignedYccTreeFN':      '%s/y.tree.aligned.ycc.%s.nwk',
            'alignedHgsnpTreeFN':    '%s/y.tree.aligned.hgSNP.%s.nwk',
            'platformYccTreeFNtp':   '%s/y.tree.platform.%%d.ycc.%s.nwk',
            'platformHgsnpTreeFNtp': '%s/y.tree.platform.%%d.hgSNP.%s.nwk',
            'bfTreeFN':              '%s/y.tree.bf.traversal.%s.txt',
            'dfTreeFN':              '%s/y.tree.df.traversal.%s.txt',
            'treeTableFN':           '%s/y.tree.table.%s.tsv',
            'bfPrimaryTreeFN':       '%s/y.tree.primary.bf.traversal.%s.txt',
            'dfPrimaryTreeFN':       '%s/y.tree.primary.df.traversal.%s.txt',
            'cleanedIsoggFN':        '%s/isogg.snps.cleaned.%s.txt',
            'uniqueIsoggFN':         '%s/isogg.snps.unique.%s.txt',
            'droppedIsoggFN':        '%s/isogg.snps.dropped.%s.txt',
        },
        'withOutdir': {
            'multiAllelicFoundFN':  '%s/multiallelic.pos',
            'pageMappingsFN':       '%s/23andMe.content.page.mappings.txt',
            'pageUpdatesFN':        '%s/23andMe.page.updates.txt',
        }
    }

    # output: haplogroup calls, log, optional files, 23andMe auxiliary files
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

    def __init__(self,
                 useDefaultCmdLineArgs=False,
                 suppressOutputAndLog=False,
                 outDir=None,
                 residList=None):
        self.residList = residList
        self.useDefaultCmdLineArgs = useDefaultCmdLineArgs
        self.suppressOutputAndLog = suppressOutputAndLog

        self.setCommandLineArgs()
        self.setParamsGeneral(outDir)
        self.validateParams()
        self.setParamsBasedOnInputType()
        self.setParamsBasedOnRunType()
        self.makeOutputDirectories()
        self.setOutputFileNamesAndOpenSome()

        if self.args.fileNamesOnly:
            self.printFileNamesAndExit()

        self.openLogAndWelcome()

        if self.runFromAblocks:
            self.setParams23andMe()
            self.get23andMeDatasets()

        if self.suppressOutputAndLog:
            self.overrideOutputGeneratingArgs()

    #----------------------------------------------------------------------
    def setParamsGeneral(self, outDir):
        'set general parameters'

        # run type: zero or one of these four will be set to True by downstream methods
        self.runFromAblocks = False
        self.runFromSampleMajorTxt = False
        self.runFromVCF = False
        self.runFromVCF4 = False
        
        # output directories
        if outDir is not None:
            self.outDir = outDir
        elif self.args.outDir is not None:
            self.outDir = self.args.outDir
        else:
            self.outDir = Config.defaultOutDir
        self.phyloOutDir = self.outDir

        # dependent args
        if self.args.example_1000Y_subset or self.args.example_single_sample_vcf:
            self.args.all_aux_output = True
        
        if self.args.all_aux_output:
            aux_output_arg_list = [
                'writeAncDerCounts', 'writeHaplogroupPaths', 'writeHaplogroupPathsDetail',
                'writeDerSNPs', 'writeDerSNPsDetail', 'writeAncSNPs', 'writeAncSNPsDetail']
            for aux_output_arg in aux_output_arg_list:
                setattr(self.args, aux_output_arg, True)
        
        # example data
        if self.args.example_1000Y_subset:
            self.args.dataFN = type(self).kgp_subset_fn
        elif self.args.example_single_sample_vcf:
            self.args.dataFN = type(self).kgp_single_sample_vcf_fn
        
        # derived 1000Y testing parameters
        self.test1000Y = (self.args.test1000Yall
                          or self.args.test1000YplatformVersion
                          or self.args.test1000Ysubset
                          or self.args.test1000YoneID)
        self.test1000YformatSpecified = (self.args.test1000Yvcf
                                         or self.args.test1000Yvcf4)

        # other parameters
        self.vcfStartCol = Config.vcfStartCol
        self.numCharsToCompare = Config.numCharsToCompareDefault
        self.prevCalledHgFN = self.args.prevCalledHgFN
        self.compareToPrevCalls = (
            self.prevCalledHgFN
            or (self.args.dataFN
                and self.args.dataFN.endswith('.resid.txt')
                and os.path.isfile(self.args.dataFN)
                and len(open(self.args.dataFN).readline().split()) > 2)
            or self.test1000Y
            or self.args.compareToMetadata)
    
    #----------------------------------------------------------------------
    def validateParams(self):
        'ensure consistency of run options'

        # preclude specification of both a data option and a test option
        if self.test1000Y:
            if self.args.dataFN:
                sys.exit('ERROR. Do not specify a 1000Y test option '
                         'and a specific data file.\n')
            if self.args.allMaleCustomers:
                sys.exit('ERROR. Do not specify a 1000Y test option '
                         'to run on all male 23andMe customers.\n')
            if self.args.primaryOnly:
                sys.exit('ERROR. Do not specify a 1000Y test option '
                         'if not importing ISOGG SNPs.\n')
        elif self.test1000YformatSpecified:
            sys.exit('ERROR. Only specify a 1000Y data format '
                     'when running on 1000Y data.\n')

        # require a resid list input file when ablock dataset specified
        if self.args.ablockDSname or self.args.ablocks_dir:
            if not self.args.dataFN or not self.args.dataFN.endswith('.resid.txt'):
                sys.exit('ERROR. Provide a resid list (-i file_label.resid.txt) '
                         'when specifying\n    a non-default ablock dataset '
                         'name or ablock directory')
            if self.args.compareToMetadata:
                sys.exit('ERROR. Will not check metadata when using '
                         'non-default ablock dataset.\n')

    #----------------------------------------------------------------------
    def setParamsBasedOnInputType(self):
        'set parameters based on input type'
        
        if self.test1000Y:
            self.setParams1000Y()
        elif self.args.dataFN:          # any type of input file, including resid list
            self.parseDataFN()
        elif self.args.allMaleCustomers:
            self.runFromAblocks = True
            self.outFNlabel = '23andMe.all.'
        elif self.residList:            # residList specified at config instantiation
            self.runFromAblocks = True
            self.outFNlabel = 'residList.'
        else:
            self.outFNlabel = ''

    def setParams1000Y(self):
        'set parameters for 1000Y testing'

        # parameters
        self.outDir = self.outDir + '.1000Y'
        self.args.alternativeRoot = 'A0-T'
        if not self.prevCalledHgFN:
            self.prevCalledHgFN = Config.thousandYhgFN

        # output-file label
        if self.args.test1000Yall:
            self.outFNlabel = '1000Y.all.'
        elif self.args.test1000Ysubset:
            self.outFNlabel = '1000Y.subset.'
        elif self.args.test1000YoneID:
            self.outFNlabel = '1000Y.%s.' % self.args.test1000YoneID
            if (self.args.singleSampleID
                and self.args.singleSampleID != self.args.test1000YoneID):
                sys.exit('ERROR. Contradiction. %s vs. %s' %
                         (self.args.singleSampleID, self.args.test1000YoneID))
        elif self.args.test1000YplatformVersion:
            self.outFNlabel = '1000Y.all.v%d.' % self.args.test1000YplatformVersion

        # input file name
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
            residFNlabel = utils.basenameNoEnding(dataFN, '.resid.txt')
            self.outFNlabel = residFNlabel
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

        if self.runFromAblocks and self.args.ablockDSname is None:
            self.outDir = self.outDir + '.23andMe'

        if self.runFromVCF4:
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

        for fn, tp in six.iteritems(Config.phyloOutputFNtpDictDict['withOutdirAndIsoggDate']):
            setattr(self, fn, tp % (self.phyloOutDir, self.isoggDate))

        for fn, tp in six.iteritems(Config.phyloOutputFNtpDictDict['withOutdir']):
            setattr(self, fn, tp % self.phyloOutDir)

        if self.args.singleSampleID:
            self.outFNlabel = '%s%s.' % (self.outFNlabel, self.args.singleSampleID)
        self.haplogroupCallsFN = self.constructOutFileName(Config.haplogroupCallsFNtp)
        self.logFN = self.constructOutFileName(Config.logFNtp)

        if self.args.writeAncDerCounts:
            self.countsAncDerFN = self.constructOutFileName(Config.countsAncDerFNtp)
        if self.args.writeHaplogroupPaths or self.args.writeHaplogroupPathsDetail:
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
            self.haplogroupRealTimeFile = open(self.haplogroupRealTimeFN, 'w', 1)

        if self.args.haplogroupToListGenotypesFor:
            self.hgGenosFN = (self.constructOutFileName(Config.hgGenosFNtp) %
                              self.args.haplogroupToListGenotypesFor)
            self.hgGenosFile = open(self.hgGenosFN, 'w', 1)

    def constructOutFileName(self, FNtp):
        'returns an output file name, given a template'

        return FNtp % (self.outDir, self.outFNlabel)

    #----------------------------------------------------------------------
    def printFileNamesAndExit(self):
        'prints input and output file names to stdout, then exits'

        print('in:  %s' % self.args.dataFN)
        print('out: %s' % self.haplogroupCallsFN)
        sys.exit()

    #----------------------------------------------------------------------
    def openLogAndWelcome(self):
        'opens log file and emits a welcome message'

        if self.suppressOutputAndLog:
            self.logFile = None
        else:
            self.logFile = open(self.logFN, 'w', 1)

        self.errAndLog('\n%s   yHaplo %s | Y-chromosome haplogroup caller\n' %
                       (utils.DASHES, __version__))
        if not self.useDefaultCmdLineArgs:
            module_name = sys.argv[0]
            if 'ttam/' in module_name:
                module_name = module_name[module_name.rfind('ttam/'):]
            else:
                module_name = module_name[module_name.rfind('yhaplo/'):]
            module_name = module_name.replace('/', '.').replace('.py', '')
            self.errAndLog('      Command: python -m %s %s\n' %
                           (module_name, ' '.join(sys.argv[1:])))
        if not self.suppressOutputAndLog:
            self.errAndLog('      Log:     %s\n' % self.logFN)
        self.errAndLog('%s' % utils.DASHES)

        self.emitWarnings()

    def errAndLog(self, message):
        'output a message to stderr and write to the log file'

        message = message.replace(Config.softwareDir + '/', '')
        sys.stderr.write(message)
        if self.logFile:
            self.logFile.write(message)

    def emitWarnings(self):
        'emit warnings for deprecated options, etc.'

        if self.args.compareToMetadata:
            self.errAndLog('\nWARNING. Deprecated option: -mdh, --compareToMetadata.\n' +
                           'The old algorithm will soon be retired.\n\n')

    #----------------------------------------------------------------------
    def setParams23andMe(self):
        'set arguments for 23andMe data'
        
        self.ablockDSname = None
        self.ablocks_dir = None
        if self.args.ablocks_dir:
            self.ablocks_dir = self.args.ablocks_dir
        elif self.args.ablockDSname:
            self.ablockDSname = self.args.ablockDSname
        else:
            self.ablockDSname = type(self).ablockDSnameDefault
        
        self.noAblocksFN = self.constructOutFileName(Config.noAblocksFNtp)
        self.noGenotypesFN = self.constructOutFileName(Config.noGenotypesFNtp)
        self.args.writeContentMappings = True
        self.numCharsToCompare = 1

    def get23andMeDatasets(self):
        '''
        get 3 datasets:
        - snp metadata
        - customer metadata (if pulling ablocks from ablockDSnameDefault)
        - ablock dataset (unless self.args.ablocks_dir has been specified)
        '''
        
        # imports
        self.errAndLog('\n%sAccessing 23andMe data...\n\n' % utils.DASHES)
        try:
            from rtk23.config import init_rtk23
            from rtk23.dataset import dataset_factory, UnknownDatasetException
            from rtk23.lib.coregen import VALUE_TO_CALL as ablockCodeToGenotypeDict
        except ImportError:
            sys.exit('ERROR. Cannot import from rtk23.\n'
                     'This run configuration is only supported '
                     'in the 23andMe research environment.')
        self.ablockCodeToGenotypeDict = ablockCodeToGenotypeDict

        # initialization
        self.errAndLog('    Initializing rtk23 ... ')
        init_rtk23()
        self.errAndLog('Done.\n')

        # SNP metadata
        self.errAndLog('    SNPs: inferring ablock indexes ... ')
        self.snpMetaDS = dataset_factory.get_dataset(Config.snpMetaDSname)
        self.setPos2ablockIndexListDict()
        self.errAndLog('Done.\n')

        # customer metadata
        if not (self.args.ablockDSname or self.ablocks_dir):
            self.errAndLog('    Customers: getting metadata...\n')
            self.customerMetaDS = dataset_factory.get_dataset(Config.customerMetaDSname)

        # genotypes
        if not self.ablocks_dir:
            self.errAndLog('    Genotypes: getting ablock dataset...\n')
            try:
                self.ablockDS = dataset_factory.get_dataset(self.ablockDSname)
            except UnknownDatasetException:
                sys.exit('\nERROR. Unknown ablock dataset: %s' % self.ablockDSname)

    def setPos2ablockIndexListDict(self):
        '''
        builds a dictionary with:
        key:   physical coordinate
        value: list of ablock indexes
        '''

        snpMetaArrayDict = self.snpMetaDS.load(Config.snpMetaColList)
        self.pos2ablockIndexListDict = defaultdict(list)
        chrom_pos_tuple_list = zip(*[snpMetaArrayDict[column]
                                     for column in Config.snpMetaColList])
        for ablockIndex, (chromosome, position) in enumerate(chrom_pos_tuple_list):
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
        self.args.writeHaplogroupPathsDetail = False
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
            self.errAndLog(
                ('Wrote genotypes at SNPs associated haplogroup %s:\n' +
                 '    %s\n\n') % (self.args.haplogroupToListGenotypesFor, self.hgGenosFN))

        if self.logFile:
            self.logFile.close()

    #----------------------------------------------------------------------
    def setCommandLineArgs(self):
        'reads command-line arguments or sets defaults if self.useDefaultCmdLineArgs'

        parser = argparse.ArgumentParser(description=DESCRIPTION,
            formatter_class=utils.RawTextWithDefaultsHelpFormatter)

        group = parser.add_argument_group('input')
        group.add_argument('-i', '--input',
            dest='dataFN', metavar='fileName',
            help='input file in one of the following formats:\n'
                 '    .resid.txt : col 1: 23andMe research IDs\n' +
                 '                 col 2: (optional) comma-separated int platforms\n' +
                 '                 col 3: (optional) previous haplogroup calls\n'
                 '    .genos.txt : sample-major text data\n'
                 '                 row 1: coordinates, col 1: sample IDs\n'
                 '    .vcf, .vcf.gz, .vcf4 : snp-major text data')

        group = parser.add_argument_group('[23andMe only] input')
        group.add_argument('-a', '--allMaleCustomers',
            dest='allMaleCustomers', action='store_true',
            help='run on all male 23andMe customers (or one, when used with -s)')
        group.add_argument('-ab', '--ablockDSname',
            dest='ablockDSname', metavar='dsName',
            help='non-default ablock dataset name')
        group.add_argument('-ad', '--ablocks_dir',
            help='directory containing ablocks,\n' +
                 'each named <resid>.ablock.npy.gz')
        
        group = parser.add_argument_group('output')
        group.add_argument('-o', '--outDir',
            dest='outDir', metavar='dirName',
            help='output directory')
       
        group = parser.add_argument_group('run on example data')
        group.add_argument('-ex', '--example_1000Y_subset',
            action='store_true',
            help='run yhaplo on a subset of 1000 Genomes data\n'
                 'and produce all auxiliary output')
        group.add_argument('-ex1', '--example_single_sample_vcf',
            action='store_true',
            help='run yhaplo on a single-sample 1000 Genomes VCF\n'
                 'and produce all auxiliary output')
        
        group = parser.add_argument_group('generate auxiliary output',
            'generate files detailing haplogroup calling for each individual')
        group.add_argument('-aao', '--all_aux_output',
            action='store_true',
            help='generate all auxilary output.\n'
                 'equivalent to these seven options:\n'
                 '--ancDerCounts --haplogroupPaths --haplogroupPathsDetail\n'
                 '--derSNPs --derSNPsDetail --ancSNPs --ancSNPsDetail')
        group.add_argument('-c', '--ancDerCounts',
            dest='writeAncDerCounts', action='store_true',
            help='counts of ancestral and derived alleles encountered\n'
                 'at each node visited (omits nodes with zero of each)')
        group.add_argument('-hp', '--haplogroupPaths',
            dest='writeHaplogroupPaths', action='store_true',
            help='sequence of branch labels from root to call,\n'
                 'with counts of derived SNPs observed')
        group.add_argument('-hpd', '--haplogroupPathsDetail',
            dest='writeHaplogroupPathsDetail', action='store_true',
            help='sequence of branch labels from root to call,\n'
                 'with counts of derived SNPs observed and lists thereof')
        group.add_argument('-ds', '--derSNPs',
            dest='writeDerSNPs', action='store_true',
            help='lists of derived SNPs on path')
        group.add_argument('-dsd', '--derSNPsDetail',
            dest='writeDerSNPsDetail', action='store_true',
            help='detailed information about each derived SNP on path')
        group.add_argument('-as', '--ancSNPs',
            dest='writeAncSNPs', action='store_true',
            help='lists of ancestral SNPs encountered in search')
        group.add_argument('-asd', '--ancSNPsDetail',
            dest='writeAncSNPsDetail', action='store_true',
            help='detailed information about each ancestral SNP\n'
                 'encountered in search')
        
        group = parser.add_argument_group('generate real-time auxilary output',
            'write haplogroup calling information as each individual is processed')
        group.add_argument('-rt', '--writeRealTime',
            dest='writeHaplogroupsRealTime', action='store_true',
            help='write haplogroups in real time. includes DFS rank,\n'
                 'to sort ex post facto: sort -nk5')
        group.add_argument('-hg', '--hgGenos',
            dest='haplogroupToListGenotypesFor', metavar='haplogroup',
            help='write genotypes observed for SNPs associated with\n'
                 'a specified node of the tree, when it is visited')
        
        group = parser.add_argument_group('traverse tree')
        group.add_argument('-b', '--breadthFirst',
            dest='traverseBF', action='store_true',
            help='write bread-first traversal')
        group.add_argument('-d', '--depthFirst',
            dest='traverseDF', action='store_true',
            help='write depth-first (pre-order) traversal')
        group.add_argument('-dt', '--depthFirstTable',
            dest='writeTreeTable', action='store_true',
            help='write depth-first (pre-order) traversal table')
        group.add_argument('-m', '--mrca', nargs=2,
            dest='mrcaHaplogroupList', metavar=('haplogroup1', 'haplogroup2'),
            help='output mrca of two haplogroups')
        group.add_argument('-sq', '--snpQuery',
            dest='querySNPname', metavar='snpName',
            help='list phylogenetic path for a query SNP')
        group.add_argument('-cm', '--contentMapping',
            dest='writeContentMappings', action='store_true',
            help='23andMe: map each node to the most recent ancestor\n'
                 'with an info page')
        group.add_argument('-pt', '--platformTrees',
            dest='writePlatformTrees', action='store_true',
            help='23andMe: write trees whose branch lengths are numbers\n'
                 'of platform sites')

        group = parser.add_argument_group(
            'compare to previously called haplogroups')
        group.add_argument('-mdh', '--compareToMetadata',
            dest='compareToMetadata', action='store_true',
            help='23andMe: compare to haplogroups called by original algorithm')
        group.add_argument('-ph', '--prevCalledHgFN',
            dest='prevCalledHgFN', metavar='fileName',
            help='import previously called haplogroups:\n'
                 'ID in first column, Haplogroup in last')
        
        group = parser.add_argument_group('change search parameters')
        group.add_argument('-ast', '--ancStopThresh',
            dest='ancStopThresh', metavar='anc_stop_thresh',
            type=int, default=ANC_STOP_THRESH_DEFAULT,
            help='BFS ancestral allele stopping condition')
        group.add_argument('-dct', '--derCollapseThresh',
            dest='derCollapseThresh', metavar='der_collapse_thresh',
            type=int, default=DER_COLLAPSE_THRESH_DEFAULT,
            help='BFS derived allele collapsing parameter')
        
        group = parser.add_argument_group('restrict input or traversal')
        group.add_argument('-po', '--primaryOnly',
            action='store_true',
            help='do NOT import ISOGG SNPs')
        group.add_argument('-r', '--root',
            dest='alternativeRoot', metavar='haplogroup',
            help='start searching tree from this branch')
        group.add_argument('-s', '--singleSample',
            dest='singleSampleID', metavar='ID',
            help='restrict to a single sample (resid for 23andMe data)')
        
        group = parser.add_argument_group('[dev only] test data')
        group.add_argument('-ta', '--test1000Yall',
            action='store_true',
            help='1000Y testing: all sites, all samples')
        group.add_argument('-ts', '--test1000Ysubset',
            action='store_true',
            help='1000Y testing: all sites, subset of samples')
        group.add_argument('-t1', '--test1000YoneID',
            metavar='ID',
            help='1000Y testing: all sites, one sample')
        group.add_argument('-tv', '--test1000YplatformSites', type=int,
            dest='test1000YplatformVersion', metavar='version',
            help='1000Y testing: 23andMe sites, all samples')
        
        group = parser.add_argument_group('[dev only] test data format')
        group.add_argument('-tvcf', '--test1000Yvcf',
            dest='test1000Yvcf', action='store_true',
            help='1000Y testing: use .vcf.gz file rather than .genos.txt')
        group.add_argument('-tvcf4', '--test1000Yvcf4',
            dest='test1000Yvcf4', action='store_true',
            help='1000Y testing: use .vcf4 file rather than .genos.txt')
        
        group = parser.add_argument_group('spot-check and exit')
        group.add_argument('-fn', '--fileNamesOnly',
            action='store_true',
            help='print file names to stdout and exit')
        group.add_argument('-v', '--version', action='version',
            version='yhaplo %s' % __version__)
        
        
        if self.useDefaultCmdLineArgs:
            self.args = parser.parse_args([])
        else:
            self.args = parser.parse_args()  # reads from sys.argv[1:]
