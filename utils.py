#!/usr/bin/python

#please edit this if needed to the location of the samtools program
samtoolsString ='samtools'

'''
Set of general utility functions for CRC Mapper
'''

#==================================================================
#===========================DEPENDENCIES===========================
#==================================================================

import os
import gzip
import time
import re
import math
from string import *

import subprocess
import datetime

from collections import defaultdict

#==================================================================
#===========================FUNCTIONS==============================
#==================================================================

'''
Locus and LocusCollection instances courtesy of Graham Ruby
'''

class Locus:
    __chrDict = dict()
    __senseDict = {'+':'+', '-':'-', '.':'.'}
    def __init__(self,chr,start,end,sense,ID='',score=0):
        coords = [start,end]
        coords.sort()
        if not(self.__chrDict.has_key(chr)): self.__chrDict[chr] = chr
        self._chr = self.__chrDict[chr]
        self._sense = self.__senseDict[sense]
        self._start = int(coords[0])
        self._end = int(coords[1])
        self._ID = ID
        self._score = score
    def ID(self): return self._ID
    def chr(self): return self._chr
    def start(self): return self._start
    def end(self): return self._end
    def len(self): return self._end - self._start + 1
    def score(self): return self._score
    def getAntisenseLocus(self):
        if self._sense=='.': return self
        else:
            switch = {'+':'-', '-':'+'}
            return Locus(self._chr,self._start,self._end,switch[self._sense])
    def coords(self): return [self._start,self._end]
    def sense(self): return self._sense
    def overlaps(self,otherLocus):
        if self.chr()!=otherLocus.chr(): return False
        elif not(self._sense=='.' or \
                 otherLocus.sense()=='.' or \
                 self.sense()==otherLocus.sense()): return False
        elif self.start() > otherLocus.end() or otherLocus.start() > self.end(): return False
        else: return True
    def contains(self,otherLocus):
        if self.chr()!=otherLocus.chr(): return False
        elif not(self._sense=='.' or \
                 otherLocus.sense()=='.' or \
                 self.sense()==otherLocus.sense()): return False
        elif self.start() > otherLocus.start() or otherLocus.end() > self.end(): return False
        else: return True
    def overlapsAntisense(self,otherLocus):
        return self.getAntisenseLocus().overlaps(otherLocus)
    def containsAntisense(self,otherLocus):
        return self.getAntisenseLocus().contains(otherLocus)
    def __hash__(self): return self._start + self._end
    def __eq__(self,other):
        if self.__class__ != other.__class__: return False
        if self.chr()!=other.chr(): return False
        if self.start()!=other.start(): return False
        if self.end()!=other.end(): return False
        if self.sense()!=other.sense(): return False
        return True
    def __ne__(self,other): return not(self.__eq__(other))
    def __str__(self): return self.chr()+'('+self.sense()+'):'+'-'.join(map(str,self.coords()))
    def plotStr(self): return self.chr() + ':' + self.sense() + ':' + '-'.join(map(str,self.coords()))
    def checkRep(self):
        pass
    def gffLine(self): return [self.chr(),self.ID(),'',self.start(),self.end(),'',self.sense(),'',self.ID()]



class LocusCollection:
    def __init__(self,loci,windowSize):
        self.__chrToCoordToLoci = dict()
        self.__loci = dict()
        self.__winSize = windowSize
        for lcs in loci: self.__addLocus(lcs)

    def __addLocus(self,lcs):
        if not(self.__loci.has_key(lcs)):
            self.__loci[lcs] = None
            if lcs.sense()=='.': chrKeyList = [lcs.chr()+'+', lcs.chr()+'-']
            else: chrKeyList = [lcs.chr()+lcs.sense()]
            for chrKey in chrKeyList:
                if not(self.__chrToCoordToLoci.has_key(chrKey)): self.__chrToCoordToLoci[chrKey] = dict()
                for n in self.__getKeyRange(lcs):
                    if not(self.__chrToCoordToLoci[chrKey].has_key(n)): self.__chrToCoordToLoci[chrKey][n] = []
                    self.__chrToCoordToLoci[chrKey][n].append(lcs)
    def __getKeyRange(self,locus):
        start = locus.start() / self.__winSize
        end = locus.end() / self.__winSize + 1
        return range(start,end)
    def __len__(self): return len(self.__loci)
        
    def append(self,new): self.__addLocus(new)
    def extend(self,newList):
        for lcs in newList: self.__addLocus(lcs)
    def hasLocus(self,locus):
        return self.__loci.has_key(locus)
    def remove(self,old):
        if not(self.__loci.has_key(old)): raise ValueError("requested locus isn't in collection")
        del self.__loci[old]
        if old.sense()=='.': senseList = ['+','-']
        else: senseList = [old.sense()]
        for k in self.__getKeyRange(old):
            for sense in senseList:
                self.__chrToCoordToLoci[old.chr()+sense][k].remove(old)
    def getWindowSize(self): return self.__winSize
    def getLoci(self): return self.__loci.keys()
    def getChrList(self):
        tempKeys = dict()
        for k in self.__chrToCoordToLoci.keys(): tempKeys[k[:-1]] = None
        return tempKeys.keys()
            
    def __subsetHelper(self,locus,sense):
        sense = sense.lower()
        if ['sense','antisense','both'].count(sense)!=1:
            raise ValueError("sense command invalid: '"+sense+"'.")
        matches = dict()
        senses = ['+','-']
        if locus.sense()=='.' or sense=='both': lamb = lambda s: True
        elif sense=='sense': lamb = lambda s: s==locus.sense()
        elif sense=='antisense': lamb = lambda s: s!=locus.sense()
        else: raise ValueError("sense value was inappropriate: '"+sense+"'.")
        for s in filter(lamb, senses):
            chrKey = locus.chr()+s
            if self.__chrToCoordToLoci.has_key(chrKey):
                for n in self.__getKeyRange(locus):
                    if self.__chrToCoordToLoci[chrKey].has_key(n):
                        for lcs in self.__chrToCoordToLoci[chrKey][n]:
                            matches[lcs] = None
        return matches.keys()
    def getOverlap(self,locus,sense='sense'):
        matches = self.__subsetHelper(locus,sense)
        realMatches = dict()
        if sense=='sense' or sense=='both':
            for i in filter(lambda lcs: lcs.overlaps(locus), matches):
                realMatches[i] = None
        if sense=='antisense' or sense=='both':
            for i in filter(lambda lcs: lcs.overlapsAntisense(locus), matches):
                realMatches[i] = None 
        return realMatches.keys()
    def getContained(self,locus,sense='sense'):
        matches = self.__subsetHelper(locus,sense)
        realMatches = dict()
        if sense=='sense' or sense=='both':
            for i in filter(lambda lcs: locus.contains(lcs), matches):
                realMatches[i] = None
        if sense=='antisense' or sense=='both':
            for i in filter(lambda lcs: locus.containsAntisense(lcs), matches):
                realMatches[i] = None
        return realMatches.keys()
    def getContainers(self,locus,sense='sense'):
        matches = self.__subsetHelper(locus,sense)
        realMatches = dict()
        if sense=='sense' or sense=='both':
            for i in filter(lambda lcs: lcs.contains(locus), matches):
                realMatches[i] = None
        if sense=='antisense' or sense=='both':
            for i in filter(lambda lcs: lcs.containsAntisense(locus), matches):
                realMatches[i] = None
        return realMatches.keys()
    def stitchCollection(self,stitchWindow=1,sense='both'):
        locusList = self.getLoci()
        oldCollection = LocusCollection(locusList,500)
        stitchedCollection = LocusCollection([],500)
        for locus in locusList:
            if oldCollection.hasLocus(locus):
                oldCollection.remove(locus)
                overlappingLoci = oldCollection.getOverlap(Locus(locus.chr(),locus.start()-stitchWindow,locus.end()+stitchWindow,locus.sense(),locus.ID()),sense)
                stitchTicker = 1
                while len(overlappingLoci) > 0:
                    stitchTicker+=len(overlappingLoci)
                    overlapCoords = locus.coords()
                    for overlappingLocus in overlappingLoci:
                        overlapCoords+=overlappingLocus.coords()
                        oldCollection.remove(overlappingLocus)
                    if sense == 'both':
                        locus = Locus(locus.chr(),min(overlapCoords),max(overlapCoords),'.',locus.ID())
                    else:
                        locus = Locus(locus.chr(),min(overlapCoords),max(overlapCoords),locus.sense(),locus.ID())
                    overlappingLoci = oldCollection.getOverlap(Locus(locus.chr(),locus.start()-stitchWindow,locus.end()+stitchWindow,locus.sense()),sense)
                locus._ID = '%s_%s_lociStitched' % (stitchTicker,locus.ID())
                stitchedCollection.append(locus)
            else:
                continue
        return stitchedCollection
    def getLoci(self): return self.__loci.keys()


'''uniquify function by Peter Bengtsson Used under a creative commons license
sourced from  here: http://www.peterbe.com/plog/uniqifiers-benchmark
'''
def uniquify(seq, idfun=None):
    if idfun is None:
        def idfun(x): return x
    seen = {}
    result = []
    for item in seq:
        marker = idfun(item)
        if marker in seen: continue
        seen[marker] = 1
        result.append(item)
    return result

'''
The under functions are courtesy of Charles Lin

'''

def parseTable(fn, sep, header = False,excel = False):
    '''takes in a table where columns are separated by a given symbol and outputs
    a nested list such that list[row][col]
    example call:
    table = parseTable('file.txt','\t')
    '''
    fh = open(fn)
    lines = fh.readlines()
    fh.close()
    if excel:
        lines = lines[0].split('\r')
    if lines[0].count('\r') > 0:
        lines = lines[0].split('\r')
    table = []
    if header == True:
        lines =lines[1:]
    for i in lines:
        table.append(i[:-1].split(sep))

    return table


def unParseTable(table, output, sep):
    '''takes in a table generated by parseTable and writes it to an output file
    takes as parameters (table, output, sep), where sep is how the file is delimited
    example call unParseTable(table, 'table.txt', '\t') for a tab del file
    '''
    fh_out = open(output,'w')
    if len(sep) == 0:
        for i in table:
            fh_out.write(str(i))
            fh_out.write('\n')
    else:
        for line in table:
            line = [str(x) for x in line]
            line = join(line,sep)

            fh_out.write(line)
            fh_out.write('\n')

    fh_out.close()

def importRefseq(refseqFile, returnMultiples = False):

    '''
    opens up a refseq file downloaded by UCSC
    '''
    refseqTable = parseTable(refseqFile,'\t')
    refseqDict = {}
    ticker = 1
    for line in refseqTable[1:]:
        if refseqDict.has_key(line[1]):
            refseqDict[line[1]].append(ticker)
        else:
            refseqDict[line[1]] = [ticker]
        ticker = ticker + 1

    multiples = []
    for i in refseqDict:
        if len(refseqDict[i]) > 1:
            multiples.append(i)

    if returnMultiples == True:
        return refseqTable,refseqDict,multiples
    else:
        return refseqTable,refseqDict

def getTSSs(geneList,refseqTable,refseqDict):
    if len(geneList) == 0:
        refseq = refseqTable
    else:
        refseq = refseqFromKey(geneList,refseqDict,refseqTable)
    TSS = []
    for line in refseq:
        if line[3] == '+':
            TSS.append(line[4])
        if line[3] == '-':
            TSS.append(line[5])
    TSS = map(int,TSS)
    return TSS

def refseqFromKey(refseqKeyList,refseqDict,refseqTable):
    typeRefseq = []
    for name in refseqKeyList:
        if refseqDict.has_key(name):
            typeRefseq.append(refseqTable[refseqDict[name][0]])
    return typeRefseq

def makeStartDict(annotFile,geneList = []):
    '''makes a dictionary keyed by refseq ID that contains information about 
    chrom/start/stop/strand/common name
    '''

    if type(geneList) == str:
        geneList = parseTable(geneList,'\t')
        geneList = [line[0] for line in geneList]
            
    if upper(annotFile).count('REFSEQ') == 1:
        refseqTable,refseqDict = importRefseq(annotFile)
        if len(geneList) == 0:
            geneList = refseqDict.keys()
        startDict = {}
        for gene in geneList:
            if refseqDict.has_key(gene) == False:
                continue
            startDict[gene]={}
            startDict[gene]['sense'] = refseqTable[refseqDict[gene][0]][3]
            startDict[gene]['chr'] = refseqTable[refseqDict[gene][0]][2]
            startDict[gene]['start'] = getTSSs([gene],refseqTable,refseqDict)
            if startDict[gene]['sense'] == '+':
                startDict[gene]['end'] =[int(refseqTable[refseqDict[gene][0]][5])]
            else:
                startDict[gene]['end'] = [int(refseqTable[refseqDict[gene][0]][4])]
            startDict[gene]['name'] = refseqTable[refseqDict[gene][0]][12]
    return startDict

def makeTSSLocus(gene,startDict,upstream,downstream):
    '''given a startDict, make a locus for any gene's TSS w/ upstream and downstream windows
    '''
    
    start = startDict[gene]['start'][0]
    if startDict[gene]['sense'] =='-':
        return Locus(startDict[gene]['chr'],start-downstream,start+upstream,'-',gene)
    else:
        return Locus(startDict[gene]['chr'],start-upstream,start+downstream,'+',gene)


def locusCollectionToGFF(locusCollection):
    lociList = locusCollection.getLoci()
    gff = []
    for locus in lociList:
        newLine = [locus.chr(),locus.ID(),'',locus.coords()[0],locus.coords()[1],'',locus.sense(),'',locus.ID()]
        gff.append(newLine)
    return gff


def makeSearchLocus(locus,upSearch,downSearch):
    '''takes a locus and expands it by a fixed upstream/downstream amount. spits out the new larger locus
    '''
    if locus.sense() == '-':
        searchLocus = Locus(locus.chr(),locus.start()-downSearch,locus.end()+upSearch,locus.sense(),locus.ID())
    else:
        searchLocus = Locus(locus.chr(),locus.start()-upSearch,locus.end()+downSearch,locus.sense(),locus.ID())
    return searchLocus


def fetchSeq(directory,chrom,start,end,UCSC=False,lineBreaks=True,header = True):
    '''function that fetches a sequence from a genome directory
    directory that contains individual chrom fasta files
    '''
    fn = directory + chrom + '.fa'
    fh = open(fn,'r')
    headerOffset = 0
    nStart = 0
    nEnd = 0
    if header:
        fh.seek(0)
        headerOffset = len(fh.readline())
    if lineBreaks:

        nStart = (start-1)/50
        nEnd = (end-1)/50
    if UCSC:
        fh.seek((start+nStart+headerOffset))
    else:
        fh.seek((start-1+nStart+headerOffset))
    span = ((end+nEnd-1)-(start+nStart-1))

    read = fh.read(span)
    if lineBreaks:
        read = read.replace('\n','')

    return read
    fh.close()

class Bam:
    '''A class for a sorted and indexed bam file that allows easy analysis of reads'''
    def __init__(self,bamFile):
        self._bam = bamFile



