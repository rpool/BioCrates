#! /usr/bin/env python

import sys
import scipy.stats
import statsmodels.stats.multitest

if(False):
    PVals = scipy.random.rand(10000,1)
    scipy.savetxt(fname='PVals.txt',
                  X=PVals,
                  fmt='%10.10e')
    scipy.savetxt(fname='pPVals.txt',
                  X=-scipy.log10(PVals),
                  fmt='%10.10e')

PVals = scipy.genfromtxt(fname='PVals.txt',
                         unpack=True)

BH = statsmodels.stats.multitest.multipletests(pvals=PVals,
                                               alpha=0.05,
                                               method='fdr_bh',
                                               returnsorted=True)
scipy.savetxt(fname='PythonBHCorrected.txt',
              X=BH[1],
              fmt='%10.10e')

ChiSq = scipy.stats.chi2.ppf((1.0-PVals),1)

scipy.savetxt(fname='Chi2.txt',
              X=ChiSq,
              fmt='%10.10e')

FChiSq = scipy.stats.chi2.pdf(ChiSq,1)

scipy.savetxt(fname='FChi2.txt',
              X=scipy.array([ChiSq,FChiSq,PVals]).T,
              fmt='%10.10e')

