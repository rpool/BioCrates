#! /usr/bin/env python

import os
import sys
import scipy

Cwd = os.getcwd()

HMPath         = os.path.join(Cwd,'HapMapB36R24.txt.gz')
HMDecompressed = os.path.join(Cwd,'HapMapB36R24.txt')
os.system('pigz -d -c -k '+HMPath+' > '+HMDecompressed)
FH          = open(HMDecompressed,'r')
HeaderList  = FH.readline().strip().split(',')
FH.close()

for Entry in HeaderList:
    if(Entry=='SNPID'):
        SNPIDColumn = HeaderList.index(Entry)
    if(Entry=='chr'):
        ChrColumn = HeaderList.index(Entry)
    if(Entry=='position'):
        PosColumn = HeaderList.index(Entry)
Arrays = scipy.loadtxt(fname=HMDecompressed,
                       dtype=str,
                       skiprows=1,
                       delimiter=',',
                       usecols=[SNPIDColumn,
                                ChrColumn,
                                PosColumn],
                       unpack=True)
os.remove(HMDecompressed)
SNPIDArray = Arrays[0]
ChrArray   = Arrays[1].astype(int)
PosArray   = Arrays[2].astype(int)

HMIndexDict = {}
for i in xrange(len(SNPIDArray)):
    HMIndexDict[SNPIDArray[i]] = i