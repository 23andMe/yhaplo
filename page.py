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

    def __init__(self, haplogroup, hgSNP):
        self.haplogroup = haplogroup
        self.hgSNP = hgSNP
        if hgSNP.find('-') > 0:
            self.snpLabel = hgSNP.split('-')[-1]
        else:
            self.snpLabel = None

    def __str__(self):
        return '%-15s %s' % (self.haplogroup, self.hgSNP)
