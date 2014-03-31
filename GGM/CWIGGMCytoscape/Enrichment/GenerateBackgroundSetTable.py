#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import scipy
import scipy.stats
import fileinput
import re

if('-g' in sys.argv):
    Lines = ''
    for Line in fileinput.input(['background.annotate.2013-12-13 22:03 PM']):
        Lines += Line
    Queries = Lines.split('////')[1:-1]

    QueryDict = {}
    for q in xrange(len(Queries)):
        QueryDict[q] = {}
        CurrentKey   = None
        for Line in Queries[q].split('\n'):
            LSplit     = Line.split('\t')
            if((len(LSplit[0].strip())>0) and
               (LSplit[0].strip()[-1]==':')):
                CurrentKey = LSplit[0].strip()[:-1]
            if((CurrentKey!=None) and
               (not QueryDict[q].has_key(CurrentKey))):
                QueryDict[q][CurrentKey] = []
                QueryDict[q][CurrentKey].append('\t'.join(LSplit[1:]))
            elif((CurrentKey!=None) and
                 (len(LSplit[1:])>0)):
                QueryDict[q][CurrentKey].append('\t'.join(LSplit[1:]))

    fw = open('BackgroundGeneSet.csv','w')
#    fw = sys.stdout
    fw.write('|'.join(['Query',\
                       'GeneHSA',\
                       'GeneSymbol(s)',\
                       'GeneEntrezID',\
                       'GOId',\
                       'GODB',\
                       'GODescr',\
                       'PathwayID',\
                       'PathWayDB',\
                       'PathwayDescr',\
                       'DiseaseID',\
                       'DiseaseDB',\
                       'DiseaseDescr'])+'\n')
    for q in QueryDict.iterkeys():
        PreStringList = []
        PreStringList.append(QueryDict[q]['Query'][0])
        StringList  = ['']*2
        try:
            for Entry in QueryDict[q]['Gene']:
                EntrySplit = Entry.split()
                StringList[0] = EntrySplit[0]
                if(len(EntrySplit)>0):
                    StringList[1] = ''.join(EntrySplit[1:])
        except(KeyError):
            pass
        PreStringList.extend(StringList)

        try:
            PreStringList.append(QueryDict[q]['Entrez Gene ID'][0])
        except(KeyError):
            PreStringList.append('')

        GOStringList = []
        if(QueryDict[q].has_key('GO')):
            for Entry in QueryDict[q]['GO']:
                EntrySplit     = Entry.split('\t')
                StringList     = ['']*3
                StringList[0]  = EntrySplit[-1]
                StringList[1]  = EntrySplit[1]
                StringList[-1] = EntrySplit[0].strip()
                GOStringList.append(StringList)

        PathwayStringList = []
        if(QueryDict[q].has_key('Pathway')):
            for Entry in QueryDict[q]['Pathway']:
                EntrySplit     = Entry.split('\t')
                StringList     = ['']*3
                StringList[0]  = EntrySplit[-1]
                StringList[1]  = EntrySplit[1]
                StringList[-1] = EntrySplit[0].strip()
                PathwayStringList.append(StringList)

        DiseaseStringList = []
        if(QueryDict[q].has_key('Disease')):
            for Entry in QueryDict[q]['Disease']:
                EntrySplit     = Entry.split('\t')
                StringList     = ['']*3
                StringList[0]  = EntrySplit[-1]
                StringList[1]  = EntrySplit[1]
                StringList[-1] = EntrySplit[0].strip()
                DiseaseStringList.append(StringList)

        MaxLength = max(len(GOStringList),len(PathwayStringList),len(DiseaseStringList))
        for i in xrange(MaxLength):
            if(i>=len(GOStringList)):
                GOStringList.append(['']*3)
            if(i>=len(PathwayStringList)):
                PathwayStringList.append(['']*3)
            if(i>=len(DiseaseStringList)):
                DiseaseStringList.append(['']*3)

        for i in xrange(MaxLength):
            String = []
            String.extend(PreStringList)
            String.extend(GOStringList[i])
            String.extend(PathwayStringList[i])
            String.extend(DiseaseStringList[i])
            fw.write('|'.join(String))
            fw.write('\n')

    fw.close()

if('-p' in sys.argv):
    fr        = open('BackgroundGeneSet.csv','r')
    BGSHeader = fr.readline().strip().split('|')
    fr.close()

    os.system('recode HTML BackgroundGeneSet.csv')

    BGSData = scipy.genfromtxt(fname='BackgroundGeneSet.csv',
                               dtype=str,
                               comments=None,
                               delimiter='|',
                               skip_header=1,
                               unpack=True)

    scipy.save(file='BackgroundGeneSet.npy',
               arr=BGSData)
    os.system('lbzip2 BackgroundGeneSet.csv')
    os.system('lbzip2 BackgroundGeneSet.npy')

if('-e' in sys.argv):
    os.system('lbzip2 -d BackgroundGeneSet.csv.bz2')
    fr        = open('BackgroundGeneSet.csv','r')
    BGSHeader = fr.readline().strip().split('|')
    fr.close()
    os.system('lbzip2 BackgroundGeneSet.csv')

#    os. system('lbzip2 -d BackgroundGeneSet.npy.bz2')
    BGSData = scipy.load('BackgroundGeneSet.npy')
#    os. system('lbzip2 BackgroundGeneSet.npy.bz2')
    print BGSHeader
    print len(BGSData),len(BGSData[0])
#   GO term enrichment
    GODBIndices       = scipy.where(BGSData[BGSHeader.index('GODB')]=='Gene Ontology')[0]
    GenesWithGOEntry  = scipy.unique(BGSData[BGSHeader.index('GeneEntrezID'),GODBIndices])
    NGenesWithGOEntry = len(GenesWithGOEntry)

    QueryGeneSymbols     = []
    fr                   = open(sys.argv[sys.argv.index('-e')+1],'r')
    for Line in fr:
        QueryGeneSymbols.append(Line.strip())
    fr.close()
    GeneSymbolArray,\
    Ind                  = scipy.unique(BGSData[BGSHeader.index('GeneSymbol(s)'),GODBIndices],return_index=True)
    GeneEntrezIDArray    = BGSData[BGSHeader.index('GeneEntrezID'),GODBIndices][Ind]
    QueryGeneEntrezIDs   = []
    for S in QueryGeneSymbols:
        boFoundS = False
        for j in xrange(len(GeneSymbolArray)):
            if(S in GeneSymbolArray[j].split(',')):
                QueryGeneEntrezIDs.append(GeneEntrezIDArray[j])
                boFoundS = True
                break
        if(not boFoundS):
            QueryGeneEntrezIDs.append('None')
    del GeneSymbolArray
    del Ind
    del GeneEntrezIDArray

    QueryGOIDs          = []
    NQueryGeneEntrezIDs = 0
    for i in xrange(len(QueryGeneEntrezIDs)):
        if(QueryGeneEntrezIDs=='None'):
            continue
        NQueryGeneEntrezIDs += 1
        Ind                  = scipy.where(BGSData[BGSHeader.index('GeneEntrezID')]==QueryGeneEntrezIDs[i])[0]
        QueryGOIDs.extend(BGSData[BGSHeader.index('GOId'),Ind].tolist())
    QueryGOIDs = scipy.unique(scipy.array(QueryGOIDs))

    for g in xrange(len(QueryGOIDs)):
        Indices                 = scipy.where(BGSData[BGSHeader.index('GOId')]==QueryGOIDs[g])[0]
        EntrezGenes             = BGSData[BGSHeader.index('GeneEntrezID'),Indices]
        NGenesWithThisGOID      = len(Indices)
        QueryGenesWithThisGOID  = scipy.intersect1d(EntrezGenes,QueryGeneEntrezIDs)
        NQueryGenesWithThisGOID = len(QueryGenesWithThisGOID)
        QueryGSymbols           = []
        for EID in QueryGenesWithThisGOID:
            QueryGSymbols.append(QueryGeneSymbols[QueryGeneEntrezIDs.index(EID)])
        EnrichmentFactor = (float(NQueryGenesWithThisGOID)/float(len(QueryGeneEntrezIDs)))/(float(NGenesWithThisGOID)/float(NGenesWithGOEntry))
        GODescr          = BGSData[BGSHeader.index('GODescr'),Indices[0]]
        PValue           = scipy.stats.hypergeom.sf(NQueryGenesWithThisGOID,NGenesWithGOEntry,NQueryGeneEntrezIDs,NGenesWithThisGOID)
        print g,\
              QueryGOIDs[g],\
              GODescr,\
              NGenesWithThisGOID,\
              NGenesWithGOEntry,\
              NQueryGenesWithThisGOID,\
              len(QueryGeneEntrezIDs),\
              ','.join(QueryGenesWithThisGOID),\
              ','.join(QueryGSymbols),\
              EnrichmentFactor,\
              PValue,\
              PValue*len(QueryGOIDs)


#    UniqueGOTerms = scipy.unique(BGSData[BGSHeader.index('GOId')])
#    print len(UniqueGOTerms)
#    for g in xrange(len(UniqueGOTerms)):
#        Indices     = scipy.where(BGSData[BGSHeader.index('GOId')]==UniqueGOTerms[g])[0]
#        NInQuerySet = 0
#        for Entry in QueryGeneEntrezIDs:
#            if(Entry in BGSData[BGSHeader.index('GeneEntrezID'),Indices].tolist()):
#                NInQuerySet += 1
#        if(NInQuerySet==0):
#            continue
#        print UniqueGOTerms[g],NInQuerySet,len(QueryGeneEntrezIDs),len(Indices),NGenesWithGOEntry

#    print QueryGeneSymbols,QueryGeneEntrezIDs
#    print GenesWithGOEntry,NGenesWithGOEntry



