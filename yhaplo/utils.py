# David Poznik
# 2016.01.05
# utils.py
#
# Defines utility functions and non-application-specific global constants.
#----------------------------------------------------------------------
from __future__ import absolute_import, print_function
import argparse
import csv
import errno
import gzip
import logging
import os
import re
import sys

#----------------------------------------------------------------------
# constants

DASHES = '-'*72 + '\n'

type2fmtDict = {bool: '%i', int: '%i', str: '%s', float: '%f'}


#----------------------------------------------------------------------
# utility functions

def basenameNoEnding(fn, ending):
    'returns the basename of a file and removes the supplied ending'

    return os.path.basename(fn)[:(1 - len(ending))]

def checkFileExistence(fn, fileDescription=None):
    'exits if file does not exist'

    if fileDescription:
        message = '%s file not found' % fileDescription
    else:
        message = 'File not found'
        
    if not os.path.isfile(fn):
        sys.exit('\nERROR. %s: %s\n' % (message, fn))

def closeFiles(fileList):
    'closes files from list, ignoring any that are set to None'
    
    for File in fileList:
        if File:
            File.close()

def compressWhitespace(myString):
    'replaces whitespace with a single space'

    return re.sub(r'\s+', ' ', myString)

def getCSVreader(inFN, delimiter='\t'):
    'opens a (possibly gzipped) file and creates a csv reader'

    extension = os.path.splitext(inFN)[1]
    if extension == '.gz':
        try: inFile = gzip.open(inFN, 'rt')
        except IOError:
            sys.exit('\nERROR. Could not open: %s\n' % inFN)
    else:
        try: inFile = open(inFN, 'r')
        except IOError:
            sys.exit('\nERROR. Could not open: %s\n' % inFN)

    return inFile, csv.reader(inFile, delimiter=delimiter)

def mkdirP(dirName):
    'makes a directory'

    try:
        os.makedirs(dirName)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(dirName):
            pass
        else:
            raise

def object2fmt(x):
    'returns a printf style format string appropriate to the object'

    return type2fmtDict[type(x)]

def printAndLogger(message):
    'output a message to stdout and to the logger'

    print(message)
    logging.info(message)

def printIterable(myIterable):
    'cycles through an iterable, printing each item'

    for item in myIterable:
        print(item)

def readPositionsSet(inFN, column = 0, logFunction = None):
    'reads positions from the specified column of a file and constructs a set'

    positionsSet = set()
    checkFileExistence(inFN, 'SNP positions')
    with open(inFN, 'r') as inFile:
        for line in inFile:
            pos = int(line.strip().split()[column])
            positionsSet.add(pos)

    message = '%5d unique positions read: %s\n' % (len(positionsSet), inFN)
    if logFunction:
        logFunction(message)
    else:
        sys.stderr.write(message)
    return positionsSet

def unimplementedMessage(methodName):
    'emits message and exits'

    sys.exit('\n\n! Unimplemented method: %s\nExiting.\n' % methodName)


#----------------------------------------------------------------------
# command-line arguments

class RawTextWithDefaultsHelpFormatter(argparse.RawDescriptionHelpFormatter):
    '''
    argparse help message formatter which:
    - retains help text formatting
    - adds default values to argument help
    combines argparse.RawTextHelpFormatter and argparse.ArgumentDefaultsHelpFormatter
    '''
    
    def _split_lines(self, text, _):
        return text.splitlines()
    
    def _get_help_string(self, action):
        help_message = action.help
        if '%(default)' not in action.help:
            if action.default is not argparse.SUPPRESS:
                defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    help_message += '\n(default: %(default)s)'
        return help_message
