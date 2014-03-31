#! /usr/bin/env python

import sys
import os
import re

fr         = open('MetaboliteClasses.txt','r')
MtbClasses = {}
for Line in fr:
    LSplit = Line.strip().split()
    MtbClasses[LSplit[1]] = LSplit[2]
fr.close()

fr = open('BioCratesTeloSizeResultsPValues.csv','r')
fr.readline()
PValsTelSize        = {}
PValsTelSizeBMICorr = {}
for Line in fr:
    LSplit = Line.strip().split(',')
    PValsTelSize[LSplit[0]]        = float(LSplit[1])
    PValsTelSizeBMICorr[LSplit[0]] = float(LSplit[2])
fr.close()

fr                = open('ggm_cn_sheet_pearson.csv','r')
MetabolitesRow    = fr.readline().strip().split(',')[1:]
CorrelationArrays = {}
for Entry in MetabolitesRow:
    CorrelationArrays[Entry] = []
MetabolitesColumn = []
for Line in fr:
    LSplit = Line.strip().split(',')
    MetabolitesColumn.append(LSplit[0])
    for i in range(len(MetabolitesRow)):
        Mtb  = MetabolitesRow[i]
        CorrelationArrays[Mtb].append(float(LSplit[i+1]))
fr.close()

for Entry in MetabolitesRow:
    if(not PValsTelSize.has_key(Entry)):
        PValsTelSize[Entry] = -1.0
    if(not PValsTelSizeBMICorr.has_key(Entry)):
        PValsTelSizeBMICorr[Entry] = -1.0

fw = open('pearson_network.txt','w')
fw.write('source'+','+\
         'target'+','\
         'interaction'+','+\
         'sourceclass'+','+\
         'correlation'+','+\
         'telspval'+','+\
         'telspvalbmicorr'+','+\
         'draw'
         '\n')
for i in range(len(MetabolitesRow)):
    Mtb = MetabolitesRow[i]

    MtbConvention = re.sub('[\(,\),\-,:,\ ]','.',Mtb)
    while(re.search('\.\.',MtbConvention)):
        MtbConvention = re.sub('\.\.','.',MtbConvention)

    for j in range(i,len(CorrelationArrays[Mtb])):
        fw.write(Mtb+','+\
                 MetabolitesColumn[j]+','+\
                 'pp'+',')
        String = 'UNKNOWN'
        if(MtbClasses.has_key(MtbConvention)):
            String = MtbClasses[MtbConvention]
        fw.write(String+',')
        fw.write(str(round(CorrelationArrays[Mtb][j],6))+',')
        fw.write(str(PValsTelSize[Mtb])+',')
        fw.write(str(PValsTelSizeBMICorr[Mtb])+',')
        if(Mtb==MetabolitesColumn[j]):
            fw.write('FALSE'+
                     '\n')
        elif(CorrelationArrays[Mtb][j]>0.1619):
            fw.write('TRUE'+
                     '\n')
        else:
            fw.write('FALSE'+
                     '\n')
fw.close()

fr                = open('ggm_cn_sheet_partial.csv','r')
MetabolitesRow    = fr.readline().strip().split(',')[1:]
CorrelationArrays = {}
for Entry in MetabolitesRow:
    CorrelationArrays[Entry] = []
MetabolitesColumn = []
for Line in fr:
    LSplit = Line.strip().split(',')
    MetabolitesColumn.append(LSplit[0])
    for i in range(len(MetabolitesRow)):
        Mtb  = MetabolitesRow[i]
        CorrelationArrays[Mtb].append(float(LSplit[i+1]))
fr.close()

fw = open('partial_network.txt','w')
fw.write('source'+','+\
         'target'+','\
         'interaction'+','+\
         'sourceclass'+','+\
         'correlation'+','+\
         'telspval'+','+\
         'telspvalbmicorr'+','+\
         'draw'+','+\
         '\n')
for i in range(len(MetabolitesRow)):
    Mtb = MetabolitesRow[i]

    MtbConvention = re.sub('[\(,\),\-,:,\ ]','.',Mtb)
    while(re.search('\.\.',MtbConvention)):
        MtbConvention = re.sub('\.\.','.',MtbConvention)

    for j in range(i,len(CorrelationArrays[Mtb])):
        fw.write(Mtb+','+\
                 MetabolitesColumn[j]+','+\
                 'pp'+',')
        String = 'UNKNOWN'
        if(MtbClasses.has_key(MtbConvention)):
            String = MtbClasses[MtbConvention]
        fw.write(String+',')
        fw.write(str(round(CorrelationArrays[Mtb][j],6))+',')
        fw.write(str(PValsTelSize[Mtb])+',')
        fw.write(str(PValsTelSizeBMICorr[Mtb])+',')
        if(Mtb==MetabolitesColumn[j]):
            fw.write('FALSE'+
                     '\n')
        elif(CorrelationArrays[Mtb][j]>0.1619):
            fw.write('TRUE'+
                     '\n')
        else:
            fw.write('FALSE'+
                     '\n')
fw.close()

