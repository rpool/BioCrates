#! /usr/bin/env python
import scipy
import sys
import re

fr         = open('HeinzResultsFdr1e-1.csv','r')
HeinzMtbs  = []
HeinzGenes = []
for Line in fr:
    LSplit = Line.strip().split(',')
    HeinzMtbs.append(re.sub('_fdr_0.1.err','',LSplit[0]))
    HeinzGenes.append(LSplit[5:])
fr.close()

for i in xrange(len(HeinzGenes)):
    HeinzGenes[i] = list(set(HeinzGenes[i]))
    try:
        del HeinzGenes[i][HeinzGenes[i].index('')]
    except(ValueError):
        pass
    HeinzGenes[i] = scipy.array(HeinzGenes[i])

fr                                  = open('Biocrates_Metabolites.csv','r')
MapDictHeinzMtbsToKORABiocratesMtbs = {}
fr.readline()
for Line in fr:
    LSplit = Line.strip().split(',')
    if(LSplit[-1]!=''):
        MapDictHeinzMtbsToKORABiocratesMtbs[LSplit[-1]] = LSplit[-2]
fr.close()

HeinzGeneDict = {}
for Key, Value in MapDictHeinzMtbsToKORABiocratesMtbs.iteritems():
    HeinzGeneDict[Value] = HeinzGenes[HeinzMtbs.index(Key)]


fr                 = open('Yoshiko.csv','r')
YoshikoClusterDict = {}
YoshikoClusters    = []
fr.readline()
for Line in fr:
    LSplit                        = Line.strip().split(',')
    YoshikoClusterDict[LSplit[0]] = LSplit[-1].split(';')
    YoshikoClusters.append(LSplit[0])
fr.close()

fw = open('MtbClusterOverlappingGenes.csv','w')
fw.write('ClusterIndex,OverlappingGenes,Metabolites\n')
for Cluster in YoshikoClusters:
    ClusterGeneDict = {}
    for Mtb in YoshikoClusterDict[Cluster]:
        ClusterGeneDict[Mtb] = HeinzGeneDict[Mtb]
    Cntr              = 0
    InterSectionArray = None
    CurrentArray      = None
    for Value in ClusterGeneDict.itervalues():
        if(Cntr==0):
            CurrentArray = Value
        InterSectionArray = scipy.intersect1d(CurrentArray,Value)
        Cntr             += 1
    fw.write(Cluster+','+';'.join(InterSectionArray.tolist())+','+';'.join(YoshikoClusterDict[str(Cluster)])+'\n')
fw.close()
