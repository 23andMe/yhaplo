#!/usr/bin/env python
#
# David Poznik
# 2016.7.26
# page.py
#
# Defines the Page class.
#----------------------------------------------------------------------

class Page(object):
    'A simple container for a 23andMe content page'

    stringTp = '%-10s %s'
    
    def __init__(self, yccOld, snpName):
        self.yccOld = yccOld
        self.snpName = snpName
        self.node = None

    def __str__(self):
        if not self.node:
            return Page.stringTp % ('.', '.')
        elif self.node.isRoot():
            return Page.stringTp % (self.node.haplogroup, self.node.haplogroup)
        else:
            return Page.stringTp % (self.node.hgSNP, self.node.haplogroup)

    def strFull(self):
        return '%-10s %-10s %s' % (self.yccOld, self.snpName, str(self))

    def setNode(self, node):
        self.node = node
