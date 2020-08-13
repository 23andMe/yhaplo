#!/usr/bin/env python
#
# David Poznik
# 2016.06.30
# convert_to_genos.py
#
# To run: python -m yhaplo.convert_to_genos
#----------------------------------------------------------------------
from __future__ import absolute_import, print_function
import argparse
import os
import sys

from . import utils

DESCRIPTION = '''
Converts data to .genos.txt format for yHaplogroup sofware.

Input format options:

1. .ped and .map

2. .23andMe.txt
    Column 1: SNP identifier (ignored)
    Column 2: Chromosome (row retained only if chromosome in CHROMOSOME_SET)
    Column 3: Physical coordinate (GRCh37 assumed)
    Column 4: Allele 1 (row retained only if allele 1 in ALLELE_SET)
    Column 5: Allele 2 (if present)
'''

CHROMOSOME_SET = {'24', 'Y'}
ALLELE_SET = set('ACGTDI')

fnEnding2fnTypeDict = {
    '.ped':         'ped',
    '.23andMe.txt': 'ttam',
    '.acom.txt':    'ttam',
}

outDir = 'converted'
utils.mkdirP(outDir)


#----------------------------------------------------------------------
# ped and map

def convertPed(pedFN, fnRoot, fnEnding):
    'reads a .ped and a .map and converts to .genos.txt'
    
    mapFN = pedFN.replace(fnEnding, '.map')
    outFN = '%s/%s.genos.txt' % (outDir, fnRoot)
    outFile = open(outFN, 'w')
    
    indexList = readMap(mapFN, outFile)
    processPed(pedFN, indexList, outFile)

    print('Output: %s\n' % outFN)
    outFile.close()

def readMap(mapFN, outFile):
    'reads a .map file'

    if not os.path.exists(mapFN):
        sys.exit('ERROR. Expecting map file: %s' % mapFN)
    print('Map: %s\n' % mapFN)
    
    positionList = list()
    indexList = list()
    index = 0
    with open(mapFN, 'r') as mapFile:
        for line in mapFile:
            chromosome, _, _, position = line.strip().split()
            if chromosome in CHROMOSOME_SET:
                positionList.append(position)
                indexList.append(index)
            index += 1
    
    outFile.write('ID\t%s\n' % '\t'.join(positionList))
    return indexList

def processPed(pedFN, indexList, outFile):
    'process a .ped file'
    
    diploidIndexList = [2 * i for i in indexList]
    numIndividuals, numFemale = 0, 0
    with open(pedFN, 'r') as inFile:
        for line in inFile:
            lineList = line.strip().split()
            sex = lineList[4]
            if sex == '2':
                numFemale += 1
                continue
            
            diploidGenoList = lineList[6:]
            haploidGenoList = list()
            for i in diploidIndexList:
                allele1, allele2 = diploidGenoList[i], diploidGenoList[i+1]
                if allele1 in ALLELE_SET and allele1 == allele2:
                    haploidGenoList.append(allele1)
                else:
                    haploidGenoList.append('.')
            
            numIndividuals += 1
            ID = '-'.join(lineList[:2])
            outFile.write('%s\t%s\n' % (ID, '\t'.join(haploidGenoList)))

    print('%5d females ignored' % numFemale)
    print('%5d individuals written' % numIndividuals)
    print('%5d markers\n' % len(indexList))


#----------------------------------------------------------------------
# 23andMe

def convertTTAM(inFN, ID):
    'reads single-sample flat format and converts to .genos.txt'

    outFN = '%s/%s.genos.txt' % (outDir, ID)
    
    genoTupleList = list()
    numNonY, numHetOrNoCall = 0, 0
    with open(inFN, 'r') as inFile:
        for line in inFile:
            if line[0] == '#' or line[:4] == 'rsid':
                continue
            
            lineList = line.strip().split()
            numFields = len(lineList)
            chromosome, position, allele1 = lineList[1:4]
            if numFields == 5:
                allele2 = lineList[4]
            elif numFields != 4:
                sys.exit('ERROR. Encountered line with %d elements:\n%s' %
                         (numFields, line))
             
            if chromosome in CHROMOSOME_SET:
                if allele1 in ALLELE_SET and (numFields == 4 or allele1 == allele2):
                    genoTupleList.append((position, allele1))
                else:
                    numHetOrNoCall += 1
            else:
                numNonY += 1
    
    with open(outFN, 'w') as outFile:
        writeLineFromTupleList(0, genoTupleList, outFile, 'ID')
        writeLineFromTupleList(1, genoTupleList, outFile, ID)
        
    print('%6d non-Y genotypes ignored' % numNonY)
    print('%6d Y-chromosome genotypes ignored (het or no-call)' % numHetOrNoCall)
    print('%6d written\n' % len(genoTupleList))
    print('Output: %s\n' % outFN)

def writeLineFromTupleList(index, tupleList, outFile, rowHeader=''):
    'given a list of tuples, writes one line with the i-th element of each tuple'
    
    outFile.write(rowHeader)
    for myTuple in tupleList:
        outFile.write('\t%s' % myTuple[index])
    
    outFile.write('\n')
    

#----------------------------------------------------------------------
# main

def main():
    parser = argparse.ArgumentParser(description=DESCRIPTION,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('inFN', type=str, help='input file name')
    args = parser.parse_args()
    inFN = args.inFN

    if not os.path.exists(inFN):
        sys.exit('ERROR. Input file does not exist: %s' % inFN)
    print('Input: %s\n' % inFN)
    
    fnType = None
    for fnEnding in fnEnding2fnTypeDict:
        if inFN.endswith(fnEnding):
            fnType = fnEnding2fnTypeDict[fnEnding]
            fnRoot = os.path.basename(inFN).replace(fnEnding, '')
            break
    
    if fnType == 'ped':
        convertPed(inFN, fnRoot, fnEnding)
    elif fnType == 'ttam':
        convertTTAM(inFN, ID=fnRoot)
    else:
        sys.exit('ERROR. Input file must be a .ped or a .23andMe.txt ' +
                 'in the corresponding format')
    
if __name__ == '__main__':
    main()
