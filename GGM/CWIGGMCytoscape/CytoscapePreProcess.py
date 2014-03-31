#! /usr/bin/env python

import sys
import os
import re

#fr         = open('MetaboliteClasses.txt','r')
#MtbClasses = {}
#for Line in fr:
#    LSplit = Line.strip().split()
#    MtbClasses[LSplit[1]] = LSplit[2]
#fr.close()

fr          = open('Biocrates_Metabolites.csv','r')
MtbClasses  = {}
ClassHeader = []
for Entry in fr.readline().strip().split(','):
    ClassHeader.append(Entry)
    MtbClasses[Entry] = []
for Line in fr:
    LSplit = Line.strip().split(',')
    for i in xrange(len(LSplit)):
        MtbClasses[ClassHeader[i]].append(LSplit[i])
fr.close()

fr                       = open('ggm_cn_sheet_pearson.csv','r')
MetabolitesRow           = fr.readline().strip().split(',')[1:]
PearsonCorrelationArrays = {}
for Entry in MetabolitesRow:
    PearsonCorrelationArrays[Entry] = {}
for Line in fr:
    LSplit = Line.strip().split(',')
    Mtb_j = LSplit[0]
    for i in range(len(MetabolitesRow)):
        Mtb_i  = MetabolitesRow[i]
        PearsonCorrelationArrays[Mtb_i][Mtb_j] = float(LSplit[i+1])
fr.close()

AllMetabolites = []
AllMetabolites.extend(MetabolitesRow)

fr                       = open('ggm_cn_sheet_partial.csv','r')
MetabolitesRow           = fr.readline().strip().split(',')[1:]
PartialCorrelationArrays = {}
for Entry in MetabolitesRow:
    PartialCorrelationArrays[Entry] = {}
for Line in fr:
    LSplit = Line.strip().split(',')
    Mtb_j  = LSplit[0]
    for i in range(len(MetabolitesRow)):
        Mtb_i  = MetabolitesRow[i]
        PartialCorrelationArrays[Mtb_i][Mtb_j] = float(LSplit[i+1])
fr.close()

AllMetabolites.extend(MetabolitesRow)
GGMMetabolites = []
GGMMetabolites.extend(AllMetabolites)
GGMMetabolites = list(set(GGMMetabolites))

fr                   = open('Jaccard.csv','r')
MetabolitesRow       = fr.readline().strip().split(',')[1:]
JaccardIndicesArrays = {}
for Entry in MetabolitesRow:
    JaccardIndicesArrays[Entry] = {}
for Line in fr:
    LSplit = Line.strip().split(',')
    Mtb_j  = LSplit[0]
    for i in range(len(MetabolitesRow)):
        Mtb_i  = MetabolitesRow[i]
        JaccardIndicesArrays[Mtb_i][Mtb_j] = float(LSplit[i+1])
fr.close()

AllMetabolites.extend(MetabolitesRow)

fr            = open('Yoshiko.csv','r')
YoshikoHeader = fr.readline().strip().split(',')
ClusterArrays = {}
for Entry in YoshikoHeader:
    ClusterArrays[Entry] = []
for Line in fr:
    LSplit = Line.strip().split(',')
    for i in xrange(len(YoshikoHeader)):
        Key = YoshikoHeader[i]
        ClusterArrays[Key].append(LSplit[i])
fr.close()

fw = open('AllNetworks.csv','w')
fw.write('KORASource'+','+\
         'KORATarget'+','\
         'KORAInteraction'+','+\
         'SharedName'+','+\
         'BioCratesSourceClass'+','+\
         'PartialCorrelation'+','+\
         'AbsPartialCorrCutoff'+','+\
         'boPartialDraw'+','+\
         'PearsonCorrelation'+','+\
         'AbsPearsonCorrCutoff'+','+\
         'boPearsonDraw'+','+\
         'JaccardIndex'+','+\
         'JaccardIndexCutoff'+','+\
         'boJaccardDraw'+','+\
         'KORASourceInENGAGEMASet'+','+\
         'KORASourceInENGAGETotalSet'+','+\
         'KORASourceInKORAGGM'+',')
String = []
for I in xrange(len(ClusterArrays['Index'])):
    String.append('KORASourceInYoshikoCluster'+ClusterArrays['Index'][I])
fw.write(','.join(String))
fw.write('\n')
AllMetabolites = list(set(AllMetabolites))
#print len(AllMetabolites)
#for Entry in AllMetabolites:
#    print Entry
#print len(GGMMetabolites), ('PC aa C38:3' in GGMMetabolites)
for i in xrange(len(AllMetabolites)):
    Mtb_i = AllMetabolites[i]

    MtbConvention = re.sub('[\(,\),\-,:,\ ,\/]','.',Mtb_i)
    while(re.search('\.\.',MtbConvention)):
        MtbConvention = re.sub('\.\.','.',MtbConvention)

    for j in xrange(len(AllMetabolites)):
        if(i==j):
            continue
        Mtb_j = AllMetabolites[j]
        fw.write(Mtb_i+','+\
                 Mtb_j+','+\
                 'pp'+',')
        fw.write(Mtb_i+' '+\
                 '(pp) '+\
                 Mtb_j+\
                 ',')
        String = 'UNKNOWN'
        try:
            String = MtbClasses['Class'][MtbClasses['KORABiocratesName'].index(Mtb_i)]
        except(ValueError):
            String = 'UNKNOWN'
        if(MtbClasses.has_key(MtbConvention)):
            String = MtbClasses[MtbConvention]
        fw.write(String+',')
        try:
            fw.write(str(round(PartialCorrelationArrays[Mtb_i][Mtb_j],6))+',')
            fw.write('0.1619,')
            fw.write(str((abs(PartialCorrelationArrays[Mtb_i][Mtb_j])>0.1619))+',')
        except(ValueError,KeyError):
            fw.write('NA,')
            fw.write('0.1619,')
            fw.write('NA,')
        try:
            fw.write(str(round(PearsonCorrelationArrays[Mtb_i][Mtb_j],6))+',')
            fw.write('0.1619,')
            fw.write(str((abs(PearsonCorrelationArrays[Mtb_i][Mtb_j])>0.1619))+',')
        except(ValueError,KeyError):
            fw.write('NA,')
            fw.write('0.1619,')
            fw.write('NA,')
        try:
            fw.write(str(round(JaccardIndicesArrays[Mtb_i][Mtb_j],6))+',')
            fw.write('0.1,')
            fw.write(str((abs(JaccardIndicesArrays[Mtb_i][Mtb_j])>0.1))+',')
        except(ValueError,KeyError):
            fw.write('NA,')
            fw.write('0.1,')
            fw.write('NA,')
        fw.write(str(JaccardIndicesArrays.has_key(Mtb_i) and \
                     PearsonCorrelationArrays.has_key(Mtb_i) and \
                     PartialCorrelationArrays.has_key(Mtb_i))+',')
        fw.write(str(Mtb_i in  MtbClasses['KORABiocratesName'])+',')
        fw.write(str(Mtb_i in  GGMMetabolites)+',')
        print Mtb_i, (Mtb_i in  GGMMetabolites)
        String = []
        for I in xrange(len(ClusterArrays['Index'])):
            String.append(str(Mtb_i in ClusterArrays['ClusterMtbs'][I].split(';')))
        fw.write(','.join(String))
        fw.write('\n')
fw.close()

