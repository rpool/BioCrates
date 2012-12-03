#! /usr/bin/env python
# python modules
import os
import argparse
import re
import scipy
import scipy.stats
import scipy.optimize
import copy
import pylab
import sys

# Homebrew modules
import Logger
import ArgumentParser
import File
import DataContainer
import Defines

#===============================================================================
# This module performs a QQ analysis on GWA output files.
#
# INPUT:
# ======
# The user should provide the following (see also "./QQ.py -h"):
# - a file that lists the metabolite names that she/he wishes to analyze
#   (option -M);
# - the full path name that contains the GWA output files (option -p);
# - the QQ analysis modes that she/he wishes to perform
#   (option -m ['P','S','PS'])
#  -> 'P' performs QQ analysis in p-value mode and stratifies on MAF and ImpQ
#  -> 'S' performs QQ analysis in chi^2 mode and stratifies on Score;
# - setting the -r flag results in a summary report in pdf format QQ analysis
#   per metabolite.
#
# OUTPUT:
# =======
# All output will be stored per metabolite and the file base name before the
# extension corresponds to the metabolite name provided in the input file by
# the user.
# The QQ plots will be stored in the 'Plots' subdirectory of the current working
# directory (CWD/Plots) in .png format.
# The stratification summaries are stored as .txt files in the CWD/Summaries
# directory.
# For option -r, the LaTeX source files and output files are stored in
# CWD/LaTeX. The directory CWD/Pdf contains the links to the generated .pdf
# summary files.
#===============================================================================

def LinearModel(X=scipy.array,
                Parameters=[]):
#   The linear regression function with ofsset 0.
    B = Parameters[0]
    Y = B*X
    return Y

def Residuals(Parameters=[],
              X=scipy.array,
              Y=scipy.array):
#   The residuals function to be called by scipy.optimize.leastsq().
    Err = Y - LinearModel(X,
                          Parameters)
    return Err

def GetPPointsArray(ArraySize=int):
    #===========================================================================
    # Copied fromm R-ppoints
    #===========================================================================
    PPointsArray   = [] # probability points array
    ThreeOverEight = 3.0/8.0
    Half           = 0.5
    N              = None
    Factor         = None
    if(ArraySize>1):
        N = ArraySize
    if(N<=10):
        Factor = ThreeOverEight
    else:
        Factor = Half
    OffSet       = Factor/float(N)
    PPointsArray = scipy.linspace(0.0+OffSet,1.0-OffSet,N)

    return PPointsArray

def LambdaEstimate(Array=scipy.array,
                   Filter=True):
    #===========================================================================
    # copied from R-GenABEL.lamdaest
    #===========================================================================
    Estimate   = None
    Ntp        = len(Array)
    QChi2Array = None

    if(Array.max()<=1.0):
#       Convert to quantile function of PValObsArray (df=1) if input are p-values.
        QChi2Array = scipy.stats.chi2.isf(Array,
                                          1)
    else:
        QChi2Array = Array

    if(Filter):
        FilterArray        = (QChi2Array>=1.0e-8)
        QChi2FilteredArray = scipy.compress(FilterArray,QChi2Array)
    else:
        QChi2FilteredArray = QChi2Array
    QChi2FilteredArray.sort()

    PPointsArray = GetPPointsArray(len(QChi2FilteredArray))
    PPointsArray = scipy.sort(scipy.stats.chi2.ppf(PPointsArray,1)) # df=1

    FilterArray          = (PPointsArray!=0.0)
    FilterArray         *= (QChi2FilteredArray!=0.0)
    PPointsArray         = scipy.compress(FilterArray,PPointsArray)
    QChi2FilteredArray   = scipy.compress(FilterArray,QChi2FilteredArray)

#   Fit PPointsArray,QChi2FilteredArray to the linear model.
    P0           = [1.0]
    PBest        = scipy.optimize.leastsq(Residuals,
                                          P0,
                                          args=(PPointsArray,QChi2FilteredArray),
                                          full_output=1,
                                          maxfev=100)
    Estimate = None
    if(type(PBest[0])==scipy.float64):
        Estimate = PBest[0]
    else:
        Estimate     = PBest[0][0]
#   Error estimation of parameter.
    Chi2 = scipy.power(PBest[2]['fvec'],2.0).sum()
    Dof  = len(QChi2FilteredArray)-len(P0)-1
    SE   = scipy.real(scipy.sqrt(PBest[1][0,0])*scipy.sqrt(Chi2/float(Dof)))

    return Estimate,SE

def PylabGetParams():
#   Set the plotting parameters: TeX mode, two columns per page in revtex mode/
    FigwidthPt  = 246.0 # pt (from revtex \showthe\columnwidth)
    InchesPerPt = 1.0/72.27
    FigWidth    = FigwidthPt*InchesPerPt
    GoldenMean  = (scipy.sqrt(5.0)-1.0)/2.0 # Aesthetic ratio
    FigHeight   = FigWidth*GoldenMean
    FigSize     = [FigWidth,FigHeight]
    Params      = {'backend': 'pdf',
                   'patch.antialiased': True,
                   'axes.labelsize': 8,
                   'axes.linewidth': 0.5,
                   'grid.color': '0.75',
                   'grid.linewidth': 0.25,
                   'grid.linestyle': ':',
                   'axes.axisbelow': False,
                   'text.fontsize': 8,
                   'legend.fontsize': 4,
                   'xtick.labelsize': 8,
                   'ytick.labelsize': 8,
                   'text.usetex': True,
                   'figure.figsize': FigSize}
    Left   = 0.16
    Bottom = 0.16
    Width  = 0.86 - Left
    Height = 0.95 - Bottom

    return Params,\
           [Left, Bottom, Width, Height]

def PlotQQFilteredOnMAF(MtbName=str,
                        PlotPath=str,
                        SummaryPath=str,
                        LPValExpArray=scipy.array,
                        LPValObsArray=scipy.array,
                        PValObsArray=scipy.array,
                        MAFArray=scipy.array,
                        MAFLevels=[],
                        RsIdArray=scipy.array,
                        Log=Logger):
#   Perform QQ analysis in p-value mode and stratify on MAF.
    PylabParameters,\
    Rectangle         = PylabGetParams()
    Size              = 2.5
    pylab.rcParams.update(PylabParameters)
    PylabFigure = pylab.figure()
    PylabFigure.clf()
    PylabAxis = PylabFigure.add_axes(Rectangle)
    PlotName  = 'QQPValModeFilteredOnMaf_'+MtbName+'.png'
    PlotName = os.path.join(PlotPath,PlotName)
#   Plot the unfiltered QQ plot.
    PylabAxis.scatter(scipy.sort(LPValExpArray),
                      scipy.sort(LPValObsArray),
                      color=Defines.Colors[0],
                      s=Size,
                      facecolor='None',
                      label=r'${\tt '+MtbName+r'}: {\rm ~all~SNPs}$')
#   Initialize the Lamba arrays and start filling them (Lambda = genomic inflation factor)
    Lambdas   = []
    SEsLambda = []
    LambdaEst,\
    SELambdaEst = LambdaEstimate(Array=PValObsArray,Filter=False)
    Lambdas.append(LambdaEst)
    SEsLambda.append(SELambdaEst)
#   Initialize the Top20 arrays and start filling them.
    Top20RsIds  = []
    Top20PVals  = []
    Top20RsIds.append([])
    Top20PVals.append([])
    ArgSortArray = scipy.argsort(PValObsArray) # Generates an index list that sorts the list.
    for i in range(20):                        # Get top 20.
        Index = ArgSortArray[i]
        Top20RsIds[-1].append(RsIdArray[Index])
        Top20PVals[-1].append(PValObsArray[Index])

    for i in range(len(MAFLevels)):
#       Plot the MAF filtered QQ plots
        Color       = Defines.Colors[i+1]
        FilterArray = None
        LabelString = None
        if(i==0):
            FilterArray = (MAFArray < MAFLevels[i])
            LabelString = '0.0 < MAF < '+str(MAFLevels[i])
        else:
            FilterArray  = (MAFArray >= MAFLevels[i-1])
            FilterArray *= (MAFArray  < MAFLevels[i])
            LabelString  = str(MAFLevels[i-1])+r'\leq {\rm MAF} <'+str(MAFLevels[i])
        ObsFP = scipy.sort(scipy.compress(FilterArray,PValObsArray))
#       Append estimates to Lambda arrays
        LambdaEst,\
        SELambdaEst = LambdaEstimate(Array=ObsFP,Filter=False)
        Lambdas.append(LambdaEst)
        SEsLambda.append(SELambdaEst)
        ObsF = scipy.sort(scipy.compress(FilterArray,LPValObsArray))
#       Generate top 20 hit list for this stratum
        Top20RsIds.append([])
        Top20PVals.append([])
        RsIdSubArray = scipy.compress(FilterArray,RsIdArray)
        PValSubArray = scipy.compress(FilterArray,PValObsArray)
        ArgSortArray = scipy.argsort(PValSubArray) # Generates an index list that sorts the list.
        for i in range(20):                        # Get top 20.
            Index = ArgSortArray[i]
            Top20RsIds[-1].append(RsIdSubArray[Index])
            Top20PVals[-1].append(PValSubArray[Index])

        scipy.random.shuffle(LPValExpArray)                          # Ensure random sample.
        ExpF = scipy.sort(scipy.compress(FilterArray,LPValExpArray)) # Sample.
        PylabAxis.scatter(ExpF,
                          ObsF,
                          color=Color,
                          s=Size,
                          facecolor='None',
                          label=r'${\tt '+MtbName+r'}: {\rm ~'+LabelString+'}$')
    MaxLPVal    = LPValObsArray.max()
    MaxLPValExp = LPValExpArray.max()
    Max         = max(MaxLPVal,MaxLPValExp)+0.5
    PylabAxis.plot([0.0,Max],
                   [0.0,Max],
                   color='black',
                   linestyle='--',
                   linewidth=0.5)
    PylabAxis.set_ylim([0.0,Max])
    PylabAxis.set_xlim([0.0,Max])
    PylabAxis.set_xlabel(r'$-\log_{10}{(P)} {\rm ~(expected)~[-]}$')
    PylabAxis.set_ylabel(r'$-\log_{10}{(P)} {\rm ~(observed)~[-]}$')
    PylabAxis.spines['right'].set_visible(False)
    PylabAxis.spines['top'].set_visible(False)
    PylabAxis.xaxis.set_ticks_position('bottom')
    PylabAxis.yaxis.set_ticks_position('left')
    Handles,Labels = PylabAxis.get_legend_handles_labels()
    PylabAxis.legend(Handles,
                     Labels,
                     fancybox=True,
                     shadow=True,
                     loc='lower right')
    PylabAxis.grid(True)
    LogString = '  ++ Saving MAF filtered plot to \"'+PlotName+'\" ...'
    print LogString
    Log.Write(LogString+'\n')
    PylabFigure.savefig(PlotName,dpi=600)
    LogString = '  -- Done ...'
    print LogString
    Log.Write(LogString+'\n')

#   Generate Lambda and top 20 summaries.
    SummaryName = 'QQPValModeFilteredOnMaf_'+MtbName+'.summary.txt'
    SummaryName = os.path.join(SummaryPath,SummaryName)

    LogString = '  ++ Saving MAF filtering summary to \"'+SummaryName+'\" ...'
    print LogString
    Log.Write(LogString+'\n')

#   Set precision
    Precision = 2
    for SE in SEsLambda:
        if(SE<1.0):
            Precision = max(Precision,int(round(-scipy.log10(SE))))
    fw = open(SummaryName,'w')
    fw.write('## Minor Allele Frequency\n')
    fw.write('# MAF level,Lambda,SELambda\n')
    fw.write('all SNPs,'+\
             str(round(Lambdas[0],Precision))+','+\
             str(round(SEsLambda[0],Precision))+\
             '\n')
    for i in range(len(MAFLevels)):
        if(i==0):
            fw.write('0.0 < MAF < '+\
                     str(MAFLevels[i])+','+\
                     str(round(Lambdas[i+1],Precision))+','+\
                     str(round(SEsLambda[i+1],Precision))+\
                     '\n')
        else:
            fw.write(str(MAFLevels[i-1])+' <= MAF < '+\
                     str(MAFLevels[i])+','+\
                     str(round(Lambdas[i+1],Precision))+','+\
                     str(round(SEsLambda[i+1],Precision))+\
                     '\n')
    fw.write('# Top 20 hits per MAF level\n')
    fw.write('#\n')
    fw.write('# MAF level: all SNPs\n')
    fw.write('# rank,rsid,p-value\n')
    for j in range(20):
        fw.write(str(j+1)+','+str(Top20RsIds[0][j])+','+str(Top20PVals[0][j])+'\n')
    for i in range(len(MAFLevels)):
        fw.write('#\n')
        if(i==0):
            fw.write('# MAF level: 0.0 < MAF < '+\
                     str(MAFLevels[i])+'\n')
        else:
            fw.write('# MAF level: '+\
                     str(MAFLevels[i-1])+' <= MAF < '+\
                     str(MAFLevels[i])+'\n')
        fw.write('# rank,rsid,p-value\n')
        for j in range(20):
            fw.write(str(j+1)+','+str(Top20RsIds[i+1][j])+','+str(Top20PVals[i+1][j])+'\n')
    fw.close()

    LogString = '  -- Done ...'
    print LogString
    Log.Write(LogString+'\n')

    return

def PlotQQFilteredOnImpQ(MtbName=str,
                         PlotPath=str,
                         SummaryPath=str,
                         LPValExpArray=scipy.array,
                         LPValObsArray=scipy.array,
                         PValObsArray=scipy.array,
                         ImpQArray=scipy.array,
                         ImpQLevels=[],
                         RsIdArray=scipy.array,
                         Log=Logger):
#   Perform QQ analysis in p-value mode and stratify on ImpQ.
    PylabParameters,\
    Rectangle         = PylabGetParams()
    Size              = 2.5
    pylab.rcParams.update(PylabParameters)
    PylabFigure = pylab.figure()
    PylabFigure.clf()
    PylabAxis = PylabFigure.add_axes(Rectangle)
    PlotName  = 'QQPValModeFilteredOnImpQ_'+MtbName+'.png'
    PlotName = os.path.join(PlotPath,PlotName)
#   Plot the unfiltered QQ plot.
    PylabAxis.scatter(scipy.sort(LPValExpArray),
                      scipy.sort(LPValObsArray),
                      color=Defines.Colors[0],
                      s=Size,
                      facecolor='None',
                      label=r'${\tt '+MtbName+r'}: {\rm ~all~SNPs}$')
#   Initialize the Lamba arrays and start filling them (Lambda = genomic inflation factor)
    Lambdas   = []
    SEsLambda = []
    LambdaEst,\
    SELambdaEst = LambdaEstimate(Array=PValObsArray,Filter=False)
    Lambdas.append(LambdaEst)
    SEsLambda.append(SELambdaEst)
#   Initialize the Top20 arrays and start filling them.
    Top20RsIds  = []
    Top20PVals  = []
    Top20RsIds.append([])
    Top20PVals.append([])
    ArgSortArray = scipy.argsort(PValObsArray) # Generates an index list that sorts the list.
    for i in range(20):                        # Get top 20.
        Index = ArgSortArray[i]
        Top20RsIds[-1].append(RsIdArray[Index])
        Top20PVals[-1].append(PValObsArray[Index])

    for i in range(len(ImpQLevels)):
#       Plot the ImpQ filtered QQ plots
        Color       = Defines.Colors[i+1]
        FilterArray = None
        LabelString = None
        if(i==0):
            FilterArray = (ImpQArray < ImpQLevels[i])
            LabelString = '0.0 < ImpQ < '+str(ImpQLevels[i])
        else:
            FilterArray  = (ImpQArray >= ImpQLevels[i-1])
            FilterArray *= (ImpQArray  < ImpQLevels[i])
            LabelString  = str(ImpQLevels[i-1])+r'\leq {\rm ImpQ} <'+str(ImpQLevels[i])
        ObsFP = scipy.sort(scipy.compress(FilterArray,PValObsArray))
#       Append estimates to Lambda arrays
        LambdaEst,\
        SELambdaEst = LambdaEstimate(Array=ObsFP,Filter=False)
        Lambdas.append(LambdaEst)
        SEsLambda.append(SELambdaEst)
        ObsF = scipy.sort(scipy.compress(FilterArray,LPValObsArray))
#       Generate top 20 hit list for this stratum
        Top20RsIds.append([])
        Top20PVals.append([])
        RsIdSubArray = scipy.compress(FilterArray,RsIdArray)
        PValSubArray = scipy.compress(FilterArray,PValObsArray)
        ArgSortArray = scipy.argsort(PValSubArray) # Generates an index list that sorts the list.
        for i in range(20):                        # Get top 20.
            Index = ArgSortArray[i]
            Top20RsIds[-1].append(RsIdSubArray[Index])
            Top20PVals[-1].append(PValSubArray[Index])
        scipy.random.shuffle(LPValExpArray)
        ExpF  = scipy.sort(scipy.compress(FilterArray,LPValExpArray))
        PylabAxis.scatter(ExpF,
                          ObsF,
                          color=Color,
                          s=Size,
                          facecolor='None',
                          label=r'${\tt '+MtbName+r'}: {\rm ~'+LabelString+'}$')
    MaxLPVal    = LPValObsArray.max()
    MaxLPValExp = LPValExpArray.max()
    Max         = max(MaxLPVal,MaxLPValExp)+0.5
    PylabAxis.plot([0.0,Max],
                   [0.0,Max],
                   color='black',
                   linestyle='--',
                   linewidth=0.5)
    PylabAxis.set_ylim([0.0,Max])
    PylabAxis.set_xlim([0.0,Max])
    PylabAxis.set_xlabel(r'$-\log_{10}{(P)} {\rm ~(expected)~[-]}$')
    PylabAxis.set_ylabel(r'$-\log_{10}{(P)} {\rm ~(observed)~[-]}$')
    PylabAxis.spines['right'].set_visible(False)
    PylabAxis.spines['top'].set_visible(False)
    PylabAxis.xaxis.set_ticks_position('bottom')
    PylabAxis.yaxis.set_ticks_position('left')
    Handles,Labels = PylabAxis.get_legend_handles_labels()
    PylabAxis.legend(Handles,
                     Labels,
                     fancybox=True,
                     shadow=True,
                     loc='lower right')
    PylabAxis.grid(True)
    LogString = '  ++ Saving imputation quality filtered plot to \"'+PlotName+'\" ...'
    print LogString
    Log.Write(LogString+'\n')
    PylabFigure.savefig(PlotName,dpi=600)
    LogString = '  -- Done ...'
    print LogString
    Log.Write(LogString+'\n')

#   Generate Lambda and top 20 summaries.
    SummaryName = 'QQPValModeFilteredOnImpQ_'+MtbName+'.summary.txt'
    SummaryName = os.path.join(SummaryPath,SummaryName)
    LogString = '  ++ Saving imputation quality filtering summary to \"'+SummaryName+'\" ...'
    print LogString
    Log.Write(LogString+'\n')

#   Set precision
    Precision = 2
    for SE in SEsLambda:
        if(SE<1.0):
            Precision = max(Precision,int(round(-scipy.log10(SE))))
    fw = open(SummaryName,'w')
    fw.write('## Imputation Quality (ImpQ)\n')
    fw.write('# ImpQ level,Lambda,SELambda\n')
    fw.write('all SNPs,'+\
             str(round(Lambdas[0],Precision))+','+\
             str(round(SEsLambda[0],Precision))+\
             '\n')
    for i in range(len(ImpQLevels)):
        if(i==0):
            fw.write(str('0.0 < ImpQ < '+\
                     str(ImpQLevels[i])+','+\
                     str(round(Lambdas[i+1],Precision))+','+\
                     str(round(SEsLambda[i+1],Precision))+\
                     '\n'))
        else:
            fw.write(str(ImpQLevels[i-1])+' <= ImpQ < '+\
                     str(ImpQLevels[i])+','+\
                     str(round(Lambdas[i+1],Precision))+','+\
                     str(round(SEsLambda[i+1],Precision))+\
                     '\n')
    fw.write('# Top 20 hits per ImpQ level\n')
    fw.write('#\n')
    fw.write('# ImpQ level: all SNPs\n')
    fw.write('# rank,rsid,p-value\n')
    for j in range(20):
        fw.write(str(j+1)+','+str(Top20RsIds[0][j])+','+str(Top20PVals[0][j])+'\n')
    for i in range(len(ImpQLevels)):
        fw.write('#\n')
        if(i==0):
            fw.write('# ImpQ level: 0.0 < ImpQ < '+\
                     str(ImpQLevels[i])+'\n')
        else:
            fw.write('# ImpQ level: '+\
                     str(ImpQLevels[i-1])+' <= ImpQ < '+\
                     str(ImpQLevels[i])+'\n')
        fw.write('# rank,rsid,p-value\n')
        for j in range(20):
            fw.write(str(j+1)+','+str(Top20RsIds[i+1][j])+','+str(Top20PVals[i+1][j])+'\n')
    fw.close()

    LogString = '  -- Done ...'
    print LogString
    Log.Write(LogString+'\n')

    return

def PlotQQFilteredOnScore(MtbName=str,
                          PlotPath=str,
                          SummaryPath=str,
                          Chi2ExpArray=scipy.array,
                          Chi2ObsArray=scipy.array,
                          ScoreArray=scipy.array,
                          ScoreLevels=[],
                          RsIdArray=scipy.array,
                          PValArray=scipy.array,
                          Log=Logger):
    PylabParameters,\
    Rectangle         = PylabGetParams()
    Size              = 2.5
    pylab.rcParams.update(PylabParameters)
    PylabFigure = pylab.figure()
    PylabFigure.clf()
    PylabAxis = PylabFigure.add_axes(Rectangle)
    PlotName  = 'QQScoreModeFilteredOnScore_'+MtbName+'.png'
    PlotName = os.path.join(PlotPath,PlotName)
#   Plot the unfiltered QQ plot.
    PylabAxis.scatter(scipy.sort(Chi2ExpArray),
                      scipy.sort(Chi2ObsArray),
                      color=Defines.Colors[0],
                      s=Size,
                      facecolor='None',
                      label=r'${\tt '+MtbName+r'}: {\rm ~all~SNPs}$')
#   Initialize the Lamba arrays and start filling them (Lambda = genomic inflation factor)
    Lambdas   = []
    SEsLambda = []
    LambdaEst,\
    SELambdaEst = LambdaEstimate(Array=Chi2ObsArray,Filter=True)
    Lambdas.append(LambdaEst)
    SEsLambda.append(SELambdaEst)
#   Initialize the Top20 arrays and start filling them.
    Top20RsIds  = []
    Top20PVals  = []
    Top20RsIds.append([])
    Top20PVals.append([])
    ArgSortArray = scipy.argsort(PValArray) # Generates an index list that sorts the list.
    for i in range(20):                     # Get top 20.
        Index = ArgSortArray[i]
        Top20RsIds[-1].append(RsIdArray[Index])
        Top20PVals[-1].append(PValArray[Index])

    for i in range(len(ScoreLevels)):
#       Plot the Score filtered QQ plots
        Color       = Defines.Colors[i+1]
        FilterArray = None
        LabelString = None
        if(i==0):
            FilterArray = (ScoreArray < ScoreLevels[i])
            LabelString = '0.0 < Score < '+str(ScoreLevels[i])
        else:
            FilterArray  = (ScoreArray >= ScoreLevels[i-1])
            FilterArray *= (ScoreArray  < ScoreLevels[i])
            LabelString  = str(ScoreLevels[i-1])+r'\leq {\rm Score} <'+str(ScoreLevels[i])
        ObsF  = scipy.sort(scipy.compress(FilterArray,Chi2ObsArray))
#       Append estimates to Lambda arrays
        LambdaEst,\
        SELambdaEst = LambdaEstimate(Array=ObsF,Filter=True)
        Lambdas.append(LambdaEst)
        SEsLambda.append(SELambdaEst)
#       Generate top 20 hit list for this stratum
        Top20RsIds.append([])
        Top20PVals.append([])
        RsIdSubArray = scipy.compress(FilterArray,RsIdArray)
        PValSubArray = scipy.compress(FilterArray,PValArray)
        ArgSortArray = scipy.argsort(PValSubArray)  # Generates an index list that sorts the list.
        for i in range(20):                         # Get top 20.
            Index = ArgSortArray[i]
            Top20RsIds[-1].append(RsIdSubArray[Index])
            Top20PVals[-1].append(PValSubArray[Index])
        scipy.random.shuffle(Chi2ExpArray)
        ExpF  = scipy.sort(scipy.compress(FilterArray,Chi2ExpArray))
        PylabAxis.scatter(ExpF,
                          ObsF,
                          color=Color,
                          s=Size,
                          facecolor='None',
                          label=r'${\tt '+MtbName+r'}: {\rm ~'+LabelString+'}$')
    MaxLPVal    = Chi2ObsArray.max()
    MaxLPValExp = Chi2ExpArray.max()
    Max         = max(MaxLPVal,MaxLPValExp)+10.0
    PylabAxis.plot([0.0,Max],
                   [0.0,Max],
                   color='black',
                   linestyle='--',
                   linewidth=0.5)
    PylabAxis.set_ylim([0.0,Max])
    PylabAxis.set_xlim([0.0,Max])
    PylabAxis.set_xlabel(r'$\chi^2 {\rm ~(expected)~[-]}$')
    PylabAxis.set_ylabel(r'$\chi^2 {\rm ~(observed)~[-]}$')
    PylabAxis.spines['right'].set_visible(False)
    PylabAxis.spines['top'].set_visible(False)
    PylabAxis.xaxis.set_ticks_position('bottom')
    PylabAxis.yaxis.set_ticks_position('left')
    Handles,Labels = PylabAxis.get_legend_handles_labels()
    PylabAxis.legend(Handles,
                     Labels,
                     fancybox=True,
                     shadow=True,
                     loc='lower right')
    PylabAxis.grid(True)
    LogString = '  ++ Saving score filtered plot to \"'+PlotName+'\" ...'
    print LogString
    Log.Write(LogString+'\n')
    PylabFigure.savefig(PlotName,dpi=600)
    LogString = '  -- Done ...'
    print LogString
    Log.Write(LogString+'\n')

#   Generate Lambda and top 20 summaries.
    SummaryName = 'QQScoreModeFilteredOnScore_'+MtbName+'.summary.txt'
    SummaryName = os.path.join(SummaryPath,SummaryName)
    LogString = '  ++ Saving score filtering summary to \"'+SummaryName+'\" ...'
    print LogString
    Log.Write(LogString+'\n')

#   Set precision
    Precision = 2
    for SE in SEsLambda:
        if(SE<1.0):
            Precision = max(Precision,int(round(-scipy.log10(SE))))
    fw = open(SummaryName,'w')
    fw.write('## Score\n')
    fw.write('# Score level,Lambda,SELambda\n')
    fw.write('all SNPs,'+\
             str(round(Lambdas[0],Precision))+','+\
             str(round(SEsLambda[0],Precision))+\
             '\n')
    for i in range(len(ScoreLevels)):
        if(i==0):
            fw.write(str('0.0 < Score < '+\
                     str(ScoreLevels[i])+','+\
                     str(round(Lambdas[i+1],Precision))+','+\
                     str(round(SEsLambda[i+1],Precision))+\
                     '\n'))
        else:
            fw.write(str(ScoreLevels[i-1])+' <= Score < '+\
                     str(ScoreLevels[i])+','+\
                     str(round(Lambdas[i+1],Precision))+','+\
                     str(round(SEsLambda[i+1],Precision))+\
                     '\n')
    fw.write('# Top 20 hits per score level\n')
    fw.write('#\n')
    fw.write('# Score level: all SNPs\n')
    fw.write('# rank,rsid,p-value\n')
    for j in range(20):
        fw.write(str(j+1)+','+str(Top20RsIds[0][j])+','+str(Top20PVals[0][j])+'\n')
    for i in range(len(ScoreLevels)):
        fw.write('#\n')
        if(i==0):
            fw.write('# Score level: 0.0 < Score < '+\
                     str(ScoreLevels[i])+'\n')
        else:
            fw.write('# Score level: '+\
                     str(ScoreLevels[i-1])+' <= Score < '+\
                     str(ScoreLevels[i])+'\n')
        fw.write('# rank,rsid,p-value\n')
        for j in range(20):
            fw.write(str(j+1)+','+str(Top20RsIds[i+1][j])+','+str(Top20PVals[i+1][j])+'\n')
    fw.close()

    LogString = '  -- Done ...'
    print LogString
    Log.Write(LogString+'\n')

    return

def QQPlotAndSummary(DCs=DataContainer.DataContainers,
                     QQModes=argparse.Namespace,
                     MtbName=str,
                     Log=Logger):
#   Initialize
    boPMode = False
    boSMode = False
    if(re.search('P',QQModes)):
        boPMode = True
    if(re.search('S',QQModes)):
        boSMode = True
    PlotPath = os.path.join(os.getcwd(),'Plots')
    if(not os.path.isdir(PlotPath)):
        os.mkdir(PlotPath)
    SummaryPath = os.path.join(os.getcwd(),'Summaries')
    if(not os.path.isdir(SummaryPath)):
        os.mkdir(SummaryPath)

    if(boPMode or boSMode):
#       Set scipy.random.seed
        scipy.random.seed(3168679810)
#       Filter out 'NA' values
        FilterArray  = (DCs.DataContainers['beta'].GetDataArray()!='NA')
        FilterArray *= (DCs.DataContainers['SE'].GetDataArray()!='NA')
        FilterArray *= (DCs.DataContainers['AF_coded_all'].GetDataArray()!='NA')
        BetaArray    = scipy.real(scipy.compress(FilterArray,DCs.DataContainers['beta'].GetDataArray()).astype(float))
        SEArray      = scipy.real(scipy.compress(FilterArray,DCs.DataContainers['SE'].GetDataArray()).astype(float))
        if((len(BetaArray)<1) or
           (len(SEArray)<1)):
            LogString  = 'xx Error message from \"QQPlotAndSummary\":\n'
            LogString += 'xx Beta and SE arrays too small after filtering on missingness!'
            print LogString
            Log.Write(LogString+'\n')
        AFArray      = scipy.real(scipy.compress(FilterArray,DCs.DataContainers['AF_coded_all'].GetDataArray()).astype(float))
        ImpQArray    = scipy.real(scipy.compress(FilterArray,DCs.DataContainers['oevar_imp'].GetDataArray()).astype(float))
        NTotArray    = scipy.real(scipy.compress(FilterArray,DCs.DataContainers['n_total'].GetDataArray()).astype(float))
        RsIdArray    = scipy.compress(FilterArray,DCs.DataContainers['SNPID'].GetDataArray())
        Chi2ObsArray = scipy.real(scipy.power((BetaArray/SEArray),2.0)).astype(float)
        FilterArray  = (Chi2ObsArray>=0.0)
        Chi2ObsArray = scipy.compress(FilterArray,Chi2ObsArray)
        PValObsArray = scipy.stats.chi2.sf(Chi2ObsArray,\
                                           1) # df=1

#       Convert AF to MAF
        CondArray    = [AFArray<=0.5]
        ChoiceArray  = [AFArray]
        MAFArray     = scipy.select(CondArray,ChoiceArray)
        CondArray    = [AFArray>0.5]
        ChoiceArray  = [1.0-AFArray]
        MAFArray    += scipy.select(CondArray,ChoiceArray)
        del CondArray
        del ChoiceArray
        del AFArray

        if(boPMode):
            MAFLevels = copy.deepcopy(Defines.MafLevels)
            MAFLevels.sort()

            ImpQLevels = copy.deepcopy(Defines.ImpLevels)
            ImpQLevels.sort()

            LPValObsArray = -scipy.log10(PValObsArray)
#           The 1's are for df=1
            PValExpArray  = scipy.stats.chi2.sf(scipy.stats.chi2.rvs(1,\
                                                                      size=len(PValObsArray)),\
                                                1)
            LPValExpArray = -scipy.log10(PValExpArray)

            LogString = '++ Plotting QQ plots in \"p-value\" mode ...'
            print LogString
            Log.Write(LogString+'\n')

            PlotQQFilteredOnMAF(MtbName,
                                PlotPath,
                                SummaryPath,
                                LPValExpArray,
                                LPValObsArray,
                                PValObsArray,
                                MAFArray,
                                MAFLevels,
                                RsIdArray,
                                Log)

            PlotQQFilteredOnImpQ(MtbName,
                                 PlotPath,
                                 SummaryPath,
                                 LPValExpArray,
                                 LPValObsArray,
                                 PValObsArray,
                                 ImpQArray,
                                 ImpQLevels,
                                 RsIdArray,
                                 Log)

            LogString = '-- Done ...'
            print LogString
            Log.Write(LogString+'\n')
        if(boSMode):
            LogString = '++ Plotting QQ plots in \"score\" mode ...'
            print LogString
            Log.Write(LogString+'\n')

#            ScoreArray  = MAFArray*ImpQArray*NTotArray*2.0  # Aaron, could you explain  this expression??
            ScoreArray  = MAFArray*ImpQArray*2.0  # Aaron, could you explain  this expression??
#            ScoreArray /= ScoreArray.max()                 # Should this be normalized to unity?
            ScoreArray *= 100.0                            # If so: normalize to 100%.
            ScoreLevels = copy.deepcopy(Defines.ScoreLevels)
            ScoreLevels.sort()
            Chi2ExpArray = scipy.stats.chi2.rvs(1,\
                                                size=len(Chi2ObsArray)) # df=1
            PlotQQFilteredOnScore(MtbName,
                                  PlotPath,
                                  SummaryPath,
                                  Chi2ExpArray,
                                  Chi2ObsArray,
                                  ScoreArray,
                                  ScoreLevels,
                                  RsIdArray,
                                  PValObsArray,
                                  Log)

            LogString = '-- Done ...'
            print LogString
            Log.Write(LogString+'\n')

    return
def LaTeXQQPModeFilteredOnMafSection(PlotFile=str,
                                     SummaryFile=str):
#   Generate the LaTeX section of the QQ analysis summary filtered on MAF
    String  = r'\section{QQ $p$--value mode: filtered on MAF}'
    String += '\n\n'
    String += r'Below the QQ plot in $p$--value mode, filtered on '
    String += r'${\rm MAF} \in [0.0,0.5]$, is shown on the left hand side. '
    String += r'The symbol $P$ denotes the $p$--value calculated from the '
    String += r'$\chi^{2}$ distribution of $\beta/SE$ with $df=1$. '
    String += r'For the observed $p$--value $\chi^{2}$ is defined as $\chi^{2}=(\beta/SE)^{2}$. '
    String += r'For the expected $p$--values, a random sample from the $\chi^{2}$ '
    String += r'distribution was taken with the same length as the observed $\chi^{2}$ array. '
    String += r'The right hand side summarizes the estimates for the genomic '
    String += r'inflation factor $\lambda_{\rm est}$ and the associated '
    String += r'standard error ${\rm SE}(\lambda_{\rm est})$, determined at '
    String += r'each MAF level.\\'
    String += '\n'
    String += r'\newline'
    String += '\n\n'
    String += r'{\tiny'
    String += '\n'
    String += r'\begin{tabular}{cc}'
    String += '\n'
    String += r'\vcent{\includegraphics[width=7.25cm,type=png,ext=.png,read=.png]{'+re.sub('\.png','',PlotFile)
    String += r'}} &'
    String += '\n'
    String += r'\begin{tabular}{l|ll|ll|ll|ll|ll|}'
    String += '\n'

#   Get the information from the SummaryFile.
#   TODO: Should be properly functionalized!
    fr                  = open(SummaryFile,'r')
    boFoundLambdaHeader = False
    LambdaDict          = {}
    SELambdaDict        = {}
    LambdaMAFLevels     = []
    boFoundTop20Header  = False
    Top20RsIdDict       = {}
    Top20PValDict       = {}
    Top20MAFLevels      = []
    for Line in fr:
        LStrip              = Line.strip()
        if(Line[0:2]=='##'):
            continue
        if(LStrip=='# MAF level,Lambda,SELambda'):
            boFoundLambdaHeader = True
            continue
        if(Line[0]=='#' and boFoundLambdaHeader):
            boFoundLambdaHeader = False
        if(boFoundLambdaHeader):
            LSplit                  = LStrip.split(',')
            LambdaDict[LSplit[0]]   = LSplit[1]
            SELambdaDict[LSplit[0]] = LSplit[2]
            LambdaMAFLevels.append(LSplit[0])
        if(LStrip=='# Top 20 hits per MAF level'):
            boFoundTop20Header = True
            continue
        if(boFoundTop20Header and re.search('# MAF level:',LStrip)):
            MAFLevel                = LStrip.split(':')[-1].strip()
            Top20RsIdDict[MAFLevel] = []
            Top20PValDict[MAFLevel] = []
            Top20MAFLevels.append(MAFLevel)
        if(len(Top20MAFLevels)!=0):
            if(LStrip[0]=='#'):
                continue
            else:
                LSplit = LStrip.split(',')
                Top20RsIdDict[Top20MAFLevels[-1]].append(LSplit[1])
                Top20PValDict[Top20MAFLevels[-1]].append(LSplit[2])
    fr.close()

    String += r'Stratum                       & \multicolumn{2}{|c|}{$a$}        & \multicolumn{2}{|c|}{$b$}             & \multicolumn{2}{|c|}{$c$}                 & \multicolumn{2}{|c|}{$d$}                & \multicolumn{2}{|c|}{$e$}          \\'
    String += '\n'
    String += r'$\lambda_{\rm est}$ '
    for Entry in LambdaMAFLevels:
        String += r'& '
        String += r'\multicolumn{2}{|c|}{'
        String += LambdaDict[Entry]
        String += r'}'
    String += r'\\'
    String += '\n'
    String += r'${\rm SE}(\lambda_{\rm est})$ '
    for Entry in LambdaMAFLevels:
        String += r'& '
        String += r'\multicolumn{2}{|c|}{'
        String += SELambdaDict[Entry]
        String += r'}'
    String += r'\\'
    String += '\n'
    String += r'[0.1cm]'
    String += '\n'
    String += r'                              & SNPID               & ${\rm p}P$ & SNPID                 & ${\rm p}P$    & SNPID                 & ${\rm p}P$        & SNPID                 & ${\rm p}P$       & SNPID                 & ${\rm p}P$ \\'
    String += '\n'
    for j in range(20):
        for i in range(len(Top20MAFLevels)):
            Entry = Top20MAFLevels[i]
            Color = 'black'
            if(j==0 and i==0):
                String += r'\multirow{20}{*}'
                String += '\n'
                String += r'{'
                String += '\n'
                String += r'\begin{sideways}'
                String += '\n'
                String += r'Top 20 hits'
                String += '\n'
                String += r'\end{sideways}'
                String += '\n'
                String += r'}'
            if(Top20RsIdDict[Entry][j] in Top20RsIdDict[Top20MAFLevels[0]]):
                Color   = 'red'
            String += r'& '
            String += r'\color{'
            String += Color
            String += r'}'
            String += r'{\tt '
            String += Top20RsIdDict[Entry][j]
            String += r'} '
            String += r'& '
            String += str(round(-scipy.log10(float(Top20PValDict[Entry][j])),2))
        String += r'\\'
        String += '\n'
    String += r'\end{tabular}\\'
    String += '\n'
    String += r'& \multicolumn{1}{p{9cm}}{Stratum symbols: '
    ChrList = [r'$a$',r'$b$',r'$c$',r'$d$',r'$e$']
    for i in range(len(LambdaMAFLevels)):
        MAFLevel  = LambdaMAFLevels[i]
        MAFLevel  = re.sub('<=',r'$\le$',MAFLevel)
        MAFLevel  = re.sub('<',r'$<$',MAFLevel)
        Chr       = ChrList[i]
        String   += Chr
        String   += r': ('
        String   += MAFLevel
        Delimiter = '), '
        if(i==len(LambdaMAFLevels)-2):
            Delimiter = ') and '
        if(i==len(LambdaMAFLevels)-1):
            Delimiter  = '). '
        String += Delimiter
    String += r'The symbol ${\rm p}P$ denotes $-\log_{10}(P)$. Note that if a SNPID is colored red, '
    String += r'it is found in the '
    String += '\"'
    String += LambdaMAFLevels[0]
    String += '\" '
    String += r'stratum.'
    String += r'}'
    String += '\n'
    String += r'\end{tabular}'
    String += '\n'
    String += r'}'
    String += '\n\n'

    return String

def LaTeXQQPModeFilteredOnImpQSection(PlotFile=str,
                                      SummaryFile=str):
#   Generate the LaTeX section of the QQ analysis summary filtered on ImpQ
    String  = r'\section{QQ $p$--value mode: filtered on ImpQ}'
    String += '\n\n'
    String += r'Below the QQ plot in $p$--value mode, filtered on imputation '
    String += r'quality ${\rm ImpQ} \in [0.0,1.0]$, is shown on the left hand side. '
    String += r'The symbol $P$ denotes the $p$--value calculated from the '
    String += r'$\chi^{2}$ distribution $df=1$. '
    String += r'For the observed $p$--value $\chi^{2}$ is defined as $\chi^{2}=(\beta/SE)^{2}$. '
    String += r'For the expected $p$--values, a random sample from the $\chi^{2}$ '
    String += r'distribution was taken with the same length as the observed $\chi^{2}$ array. '
    String += r'The right hand side summarizes the estimates for the genomic '
    String += r'inflation factor $\lambda_{\rm est}$ and the associated '
    String += r'standard error ${\rm SE}(\lambda_{\rm est})$, determined at '
    String += r'each ImpQ level.\\'
    String += '\n'
    String += r'\newline'
    String += '\n'
    String += r'\newline'
    String += '\n\n'
    String += r'{\tiny'
    String += '\n'
    String += r'\begin{tabular}{cc}'
    String += '\n'
    String += r'\vcent{\includegraphics[width=7.25cm,type=png,ext=.png,read=.png]{'+re.sub('\.png','',PlotFile)
    String += r'}} &'
    String += '\n'
    String += r'\begin{tabular}{l|ll|ll|ll|ll|ll|}'
    String += '\n'

#   Get the information from the SummaryFile.
#   TODO: Should be properly functionalized!
    fr                  = open(SummaryFile,'r')
    boFoundLambdaHeader = False
    LambdaDict          = {}
    SELambdaDict        = {}
    LambdaImpQLevels    = []
    boFoundTop20Header  = False
    Top20RsIdDict       = {}
    Top20PValDict       = {}
    Top20ImpQLevels     = []
    for Line in fr:
        LStrip              = Line.strip()
        if(Line[0:2]=='##'):
            continue
        if(LStrip=='# ImpQ level,Lambda,SELambda'):
            boFoundLambdaHeader = True
            continue
        if(Line[0]=='#' and boFoundLambdaHeader):
            boFoundLambdaHeader = False
        if(boFoundLambdaHeader):
            LSplit                  = LStrip.split(',')
            LambdaDict[LSplit[0]]   = LSplit[1]
            SELambdaDict[LSplit[0]] = LSplit[2]
            LambdaImpQLevels.append(LSplit[0])
        if(LStrip=='# Top 20 hits per ImpQ level'):
            boFoundTop20Header = True
            continue
        if(boFoundTop20Header and re.search('# ImpQ level:',LStrip)):
            ImpQLevel                = LStrip.split(':')[-1].strip()
            Top20RsIdDict[ImpQLevel] = []
            Top20PValDict[ImpQLevel] = []
            Top20ImpQLevels.append(ImpQLevel)
        if(len(Top20ImpQLevels)!=0):
            if(LStrip[0]=='#'):
                continue
            else:
                LSplit = LStrip.split(',')
                Top20RsIdDict[Top20ImpQLevels[-1]].append(LSplit[1])
                Top20PValDict[Top20ImpQLevels[-1]].append(LSplit[2])
    fr.close()

    String += r'Stratum                       & \multicolumn{2}{|c|}{$a$}        & \multicolumn{2}{|c|}{$b$}             & \multicolumn{2}{|c|}{$c$}                 & \multicolumn{2}{|c|}{$d$}                & \multicolumn{2}{|c|}{$e$}          \\'
    String += '\n'
    String += r'$\lambda_{\rm est}$ '
    for Entry in LambdaImpQLevels:
        String += r'& '
        String += r'\multicolumn{2}{|c|}{'
        String += LambdaDict[Entry]
        String += r'}'
    String += r'\\'
    String += '\n'
    String += r'${\rm SE}(\lambda_{\rm est})$ '
    for Entry in LambdaImpQLevels:
        String += r'& '
        String += r'\multicolumn{2}{|c|}{'
        String += SELambdaDict[Entry]
        String += r'}'
    String += r'\\'
    String += '\n'
    String += r'[0.1cm]'
    String += '\n'
    String += r'                              & SNPID               & ${\rm p}P$ & SNPID                 & ${\rm p}P$    & SNPID                 & ${\rm p}P$        & SNPID                 & ${\rm p}P$       & SNPID                 & ${\rm p}P$ \\'
    String += '\n'
    for j in range(20):
        for i in range(len(Top20ImpQLevels)):
            Entry = Top20ImpQLevels[i]
            Color = 'black'
            if(j==0 and i==0):
                String += r'\multirow{20}{*}'
                String += '\n'
                String += r'{'
                String += '\n'
                String += r'\begin{sideways}'
                String += '\n'
                String += r'Top 20 hits'
                String += '\n'
                String += r'\end{sideways}'
                String += '\n'
                String += r'}'
            if(Top20RsIdDict[Entry][j] in Top20RsIdDict[Top20ImpQLevels[0]]):
                Color   = 'red'
            String += r'& '
            String += r'\color{'
            String += Color
            String += r'}'
            String += r'{\tt '
            String += Top20RsIdDict[Entry][j]
            String += r'} '
            String += r'& '
            String += str(round(-scipy.log10(float(Top20PValDict[Entry][j])),2))
        String += r'\\'
        String += '\n'
    String += r'\end{tabular}\\'
    String += '\n'
    String += r'& \multicolumn{1}{p{9cm}}{Stratum symbols: '
    ChrList = [r'$a$',r'$b$',r'$c$',r'$d$',r'$e$']
    for i in range(len(LambdaImpQLevels)):
        ImpQLevel = LambdaImpQLevels[i]
        ImpQLevel = re.sub('<=',r'$\le$',ImpQLevel)
        ImpQLevel = re.sub('<',r'$<$',ImpQLevel)
        Chr       = ChrList[i]
        String   += Chr
        String   += r': ('
        String   += ImpQLevel
        Delimiter = '), '
        if(i==len(LambdaImpQLevels)-2):
            Delimiter = ') and '
        if(i==len(LambdaImpQLevels)-1):
            Delimiter  = '). '
        String += Delimiter
    String += r'The symbol ${\rm p}P$ denotes $-\log_{10}(P)$. Note that if a SNPID is colored red, '
    String += r'it is found in the '
    String += '\"'
    String += LambdaImpQLevels[0]
    String += '\" '
    String += r'stratum.'
    String += r'}'
    String += '\n'
    String += r'\end{tabular}'
    String += '\n'
    String += r'}'
    String += '\n\n'

    return String

def LaTeXQQScoreModeFilteredQSection(PlotFile=str,
                                     SummaryFile=str):
#   Generate the LaTeX section of the QQ analysis summary filtered on Score
    String  = r'\section{QQ $\chi^{2}$ mode: filtered on score}'
    String += '\n\n'
    String += r'Below the QQ plot in $\chi^{2}$ mode, filtered on score '
    String += r'quality ${\rm Score} \in [0.0,100.0]$, is shown on the left hand side. '
    String += r'The observed $\chi^{2}$ values ($df=1$) were determined by squaring '
    String += r'$\beta/SE$. '
    String += r'The right hand side summarizes the estimates for the genomic '
    String += r'inflation factor $\lambda_{\rm est}$ and the associated '
    String += r'standard error ${\rm SE}(\lambda_{\rm est})$, determined at '
    String += r'each score level.\\'
    String += '\n\n'
    String += r'{\tiny'
    String += '\n'
    String += r'\begin{tabular}{cc}'
    String += '\n'
    String += r'\vcent{\includegraphics[width=7.25cm,type=png,ext=.png,read=.png]{'+re.sub('\.png','',PlotFile)
    String += r'}} &'
    String += '\n'
    String += r'\begin{tabular}{l|ll|ll|ll|ll|ll|}'
    String += '\n'

#   Get the information from the SummaryFile.
#   TODO: Should be properly functionalized!
    fr                  = open(SummaryFile,'r')
    boFoundLambdaHeader = False
    LambdaDict          = {}
    SELambdaDict        = {}
    LambdaScoreLevels   = []
    boFoundTop20Header  = False
    Top20RsIdDict       = {}
    Top20PValDict       = {}
    Top20ScoreLevels    = []
    for Line in fr:
        LStrip              = Line.strip()
        if(Line[0:2]=='##'):
            continue
        if(LStrip=='# Score level,Lambda,SELambda'):
            boFoundLambdaHeader = True
            continue
        if(Line[0]=='#' and boFoundLambdaHeader):
            boFoundLambdaHeader = False
        if(boFoundLambdaHeader):
            LSplit                  = LStrip.split(',')
            LambdaDict[LSplit[0]]   = LSplit[1]
            SELambdaDict[LSplit[0]] = LSplit[2]
            LambdaScoreLevels.append(LSplit[0])
        if(LStrip=='# Top 20 hits per score level'):
            boFoundTop20Header = True
            continue
        if(boFoundTop20Header and re.search('# Score level:',LStrip)):
            ScoreLevel                = LStrip.split(':')[-1].strip()
            Top20RsIdDict[ScoreLevel] = []
            Top20PValDict[ScoreLevel] = []
            Top20ScoreLevels.append(ScoreLevel)
        if(len(Top20ScoreLevels)!=0):
            if(LStrip[0]=='#'):
                continue
            else:
                LSplit = LStrip.split(',')
                Top20RsIdDict[Top20ScoreLevels[-1]].append(LSplit[1])
                Top20PValDict[Top20ScoreLevels[-1]].append(LSplit[2])
    fr.close()

    String += r'Stratum                       & \multicolumn{2}{|c|}{$a$}        & \multicolumn{2}{|c|}{$b$}             & \multicolumn{2}{|c|}{$c$}                 & \multicolumn{2}{|c|}{$d$}                & \multicolumn{2}{|c|}{$e$}          \\'
    String += '\n'
    String += r'$\lambda_{\rm est}$ '
    for Entry in LambdaScoreLevels:
        String += r'& '
        String += r'\multicolumn{2}{|c|}{'
        String += LambdaDict[Entry]
        String += r'}'
    String += r'\\'
    String += '\n'
    String += r'${\rm SE}(\lambda_{\rm est})$ '
    for Entry in LambdaScoreLevels:
        String += r'& '
        String += r'\multicolumn{2}{|c|}{'
        String += SELambdaDict[Entry]
        String += r'}'
    String += r'\\'
    String += '\n'
    String += r'[0.1cm]'
    String += '\n'
    String += r'                              & SNPID               & ${\rm p}P$ & SNPID                 & ${\rm p}P$    & SNPID                 & ${\rm p}P$        & SNPID                 & ${\rm p}P$       & SNPID                 & ${\rm p}P$ \\'
    String += '\n'
    for j in range(20):
        for i in range(len(Top20ScoreLevels)):
            Entry = Top20ScoreLevels[i]
            Color = 'black'
            if(j==0 and i==0):
                String += r'\multirow{20}{*}'
                String += '\n'
                String += r'{'
                String += '\n'
                String += r'\begin{sideways}'
                String += '\n'
                String += r'Top 20 hits'
                String += '\n'
                String += r'\end{sideways}'
                String += '\n'
                String += r'}'
            if(Top20RsIdDict[Entry][j] in Top20RsIdDict[Top20ScoreLevels[0]]):
                Color   = 'red'
            String += r'& '
            String += r'\color{'
            String += Color
            String += r'}'
            String += r'{\tt '
            String += Top20RsIdDict[Entry][j]
            String += r'} '
            String += r'& '
            String += str(round(-scipy.log10(float(Top20PValDict[Entry][j])),2))
        String += r'\\'
        String += '\n'
    String += r'\end{tabular}\\'
    String += '\n'
    String += r'& \multicolumn{1}{p{9cm}}{Stratum symbols: '
    ChrList = [r'$a$',r'$b$',r'$c$',r'$d$',r'$e$']
    for i in range(len(LambdaScoreLevels)):
        ScoreLevel = LambdaScoreLevels[i]
        ScoreLevel = re.sub('<=',r'$\le$',ScoreLevel)
        ScoreLevel = re.sub('<',r'$<$',ScoreLevel)
        Chr        = ChrList[i]
        String    += Chr
        String    += r': ('
        String    += ScoreLevel
        Delimiter  = '), '
        if(i==len(LambdaScoreLevels)-2):
            Delimiter = ') and '
        if(i==len(LambdaScoreLevels)-1):
            Delimiter  = '). '
        String += Delimiter
    String += r'The symbol ${\rm p}P$ denotes $-\log_{10}(P)$. Note that if a SNPID is colored red, '
    String += r'it is found in the '
    String += '\"'
    String += LambdaScoreLevels[0]
    String += '\" '
    String += r'stratum.'
    String += r'}'
    String += '\n'
    String += r'\end{tabular}'
    String += '\n'
    String += r'}'
    String += '\n\n'
    return String

def LaTeXPreamble():
#   Generate the LaTeX preamble
    String  =r'\documentclass[pre,amsmath,onecolumn,floatfix,fleqn,a4paper,superscriptaddress]{revtex4}'
    String += '\n'
    String +=r'\usepackage{latexsym}'
    String += '\n'
    String +=r'\usepackage{amsmath,amsfonts,amssymb}'
    String += '\n'
    String +=r'\usepackage{graphicx}'
    String += '\n'
    String +=r'\usepackage{subfigure}'
    String += '\n'
    String +=r'\usepackage{float}'
    String += '\n'
    String +=r'\usepackage{pifont}'
    String += '\n'
    String +=r'\usepackage{color}'
    String += '\n'
    String +=r'\usepackage{multirow}'
    String += '\n'
    String +=r'\usepackage{rotating}'
    String += '\n'
    String += '\n'
    String +=r'\renewcommand{\textfraction}{0.00} \renewcommand{\topfraction}{1.0}'
    String += '\n'
    String +=r'\renewcommand{\bottomfraction}{1.0}'
    String += '\n'
    String +=r'\renewcommand{\floatpagefraction}{1}'
    String += '\n'
    String +=r'\newcommand\vcent[1]{\ensuremath{\vcenter{\hbox{{#1}}}}}'
    String += '\n'
    String += '\n'
    return String

def LaTeXTitle(MtbName=str):
#   Set and return the title of the report.
    String  = r'\title{QQ Analysis Report for Metabolite '
    String += r'{\tt '
    String += MtbName
    String += r'}'
    String += r'}'
    String += '\n'
    return String

def LaTeXDate():
#   Set and return the date of the report.
    String  = r'\date{\today}'
    String += '\n'
    String += '\n'
    return String

def LaTeXBeginDocument():
#   Begin document
    String  = r'\begin{document}'
    String += '\n'
    return String

def LaTeXMakeTitle():
#   Make title
    String  = r'\maketitle'
    String += '\n'
    String += '\n'
    return String

def LaTeXEndDocument():
#   End document
    String  = r'\end{document}'
    String += '\n'
    return String


def GenerateLaTeXReport(Log=Logger,
                        MtbNames=[]):
#   Initialize and log
    LogString = '++ Generating QQ analysis reports using LaTeX ...'
    print LogString
    Log.Write(LogString+'\n')
    BasePath  = os.getcwd()
    LaTeXPath = os.path.join(BasePath,'LaTeX')
    if(not os.path.isdir(LaTeXPath)):
        os.mkdir(LaTeXPath)
    PdfPath = os.path.join(BasePath,'Pdf')
    if(not os.path.isdir(PdfPath)):
        os.mkdir(PdfPath)

#   Only output sections for which there is a plot && for which there is a summary file.
    PlotPath    = os.path.join(BasePath,'Plots')
    SummaryPath = os.path.join(BasePath,'Summaries')
    for MtbName in MtbNames:
        LaTeXSrcFile    = os.path.join(LaTeXPath,'QQReport_'+MtbName+'.tex')
        LaTeXStdoutFile = os.path.join(LaTeXPath,'QQReport_'+MtbName+'.stdout')
        LaTeXStdErrFile = os.path.join(LaTeXPath,'QQReport_'+MtbName+'.stderr')

        LogString  = '  ++ Generating LaTeX source file for metabolite \"'+MtbName+'\" ...\n'
        LogString += '     LaTeX source file: \"'+LaTeXSrcFile+'\" ...'
        print LogString
        Log.Write(LogString+'\n')
        fw = open(LaTeXSrcFile,'w')
        fw.write(LaTeXPreamble())
        fw.write(LaTeXDate())
        fw.write(LaTeXBeginDocument())
        fw.write(LaTeXTitle(MtbName))
        fw.write(LaTeXMakeTitle())
#       TODO: should be properly functionalized!
        boHavePlot       = False
        QQMafPlotFile    = None
        boHaveSummary    = False
        QQMafSummaryFile = None
        for FileName in os.listdir(PlotPath):
            if(FileName=='QQPValModeFilteredOnMaf_'+MtbName+'.png'):
                QQMafPlotFile = os.path.join(PlotPath,FileName)
                boHavePlot    = True
        for FileName in os.listdir(SummaryPath):
            if(FileName=='QQPValModeFilteredOnMaf_'+MtbName+'.summary.txt'):
                QQMafSummaryFile = os.path.join(SummaryPath,FileName)
                boHaveSummary    = True
        if(boHavePlot and boHaveSummary):
            fw.write(LaTeXQQPModeFilteredOnMafSection(QQMafPlotFile,
                                                      QQMafSummaryFile))
#       TODO: should be properly functionalized!
        boHavePlot        = False
        QQImpQPlotFile    = None
        boHaveSummary     = False
        QQImpQSummaryFile = None
        for FileName in os.listdir(PlotPath):
            if(FileName=='QQPValModeFilteredOnImpQ_'+MtbName+'.png'):
                QQImpQPlotFile = os.path.join(PlotPath,FileName)
                boHavePlot     = True
        for FileName in os.listdir(SummaryPath):
            if(FileName=='QQPValModeFilteredOnImpQ_'+MtbName+'.summary.txt'):
                QQImpQSummaryFile = os.path.join(SummaryPath,FileName)
                boHaveSummary     = True
        if(boHavePlot and boHaveSummary):
            fw.write(LaTeXQQPModeFilteredOnImpQSection(QQImpQPlotFile,
                                                       QQImpQSummaryFile))
#       TODO: should be properly functionalized!
        boHavePlot         = False
        QQScorePlotFile    = None
        boHaveSummary      = False
        QQScoreSummaryFile = None
        for FileName in os.listdir(PlotPath):
            if(FileName=='QQScoreModeFilteredOnScore_'+MtbName+'.png'):
                QQScorePlotFile = os.path.join(PlotPath,FileName)
                boHavePlot      = True
        for FileName in os.listdir(SummaryPath):
            if(FileName=='QQScoreModeFilteredOnScore_'+MtbName+'.summary.txt'):
                QQScoreSummaryFile = os.path.join(SummaryPath,FileName)
                boHaveSummary      = True
        if(boHavePlot and boHaveSummary):
            fw.write(LaTeXQQScoreModeFilteredQSection(QQScorePlotFile,
                                                      QQScoreSummaryFile))
        fw.write(LaTeXEndDocument())
        fw.close()
        LogString = '  -- Done ...'
        print LogString
        Log.Write(LogString+'\n')

        PdfOutFile = os.path.join(PdfPath,'QQReport_'+MtbName+'.pdf')

        LogString  = '  ++ Generating pdf summary source file using \"pdflatex\" for metabolite \"'+MtbName+'\" ...\n'
        LogString += '     Pdf output file: \"'+PdfOutFile+'\" ...'
        print LogString
        Log.Write(LogString+'\n')
#       The .tex source is now created, let's make a document of it :-)
        os.chdir(LaTeXPath)
        os.system('pdflatex '+LaTeXSrcFile+' > '+LaTeXStdoutFile+' 2> '+LaTeXStdErrFile)
        os.chdir(BasePath)
        os.system('ln -sf '+re.sub('.tex','.pdf',LaTeXSrcFile)+' '+PdfOutFile)
        LogString = '  -- Done ...'
        print LogString
        Log.Write(LogString+'\n')

    LogString = '-- Done ...'
    print LogString
    Log.Write(LogString+'\n')

    return

def main(ExecutableName):

    #===========================================================================
    # START Initialization
    #===========================================================================
    ArgParser,\
    Arguments   = ArgumentParser.ParseArguments()

    Ext = ''
    Log = Logger.Logger(ExecutableName,
                        Ext)
    LogString  = '## START TIMESTAMP\n'
    LogString += str(Log.GetStartDate())+'\n'
    LogString += '## END TIMESTAMP'
    print LogString
    Log.Write(LogString+'\n')
    LogString = Log.GetStartLogString()
    print LogString
    Log.Write(LogString+'\n')
    ArgumentParser.LogArguments(Log,
                                ArgParser,
                                Arguments)
    #===========================================================================
    # END Initialization
    #===========================================================================

    #===========================================================================
    # START Do the work!
    #===========================================================================
    # Parse Arguments.GwaFiles
    LogString = '++ Parsing \"'+Arguments.MtbNameFile+'\" ...'
    print LogString
    Log.Write(LogString+'\n')
    MtbNameFile = File.File(Name=Arguments.MtbNameFile,
                            boHeader=False)
    MtbNameFile.SetFileHandle(Mode='r')
    MtbNameFileDCs = MtbNameFile.ParseToDataContainers()
    MtbNameFile.Close()
    MtbNameFile.Cleanup()
    del MtbNameFile
    LogString = '-- Done ...'
    print LogString
    Log.Write(LogString+'\n')

    # Loop over GwaFiles
    MtbNames        = []
    GwaDataFileList = os.listdir(Arguments.GWADataPath)
    for MtbName in MtbNameFileDCs.DataContainers['0'].GetDataArray(): # '0' because there should be ONE column && NO header
        MtbNames.append(MtbName)
        GwaFileName = None
        DelIndex    = None
        for FileName in GwaDataFileList:
            FSplit      = FileName.split('_')
            GwaFileName = None
            if(len(FSplit)>3):
                Name = FSplit[2]
            if(Name==MtbName):
                DelIndex    = GwaDataFileList.index(FileName)
                GwaFileName = os.path.join(Arguments.GWADataPath,FileName)
                break
        del GwaDataFileList[DelIndex] # reduce size of GwaDataFileList
#       Parse the GWA output file that corresponds to the MtbName
        LogString = '++ Parsing \"'+GwaFileName+'\" ...'
        print LogString
        Log.Write(LogString+'\n')
        GwaFile = File.File(Name=GwaFileName,
                            boHeader=True)
        GwaFile.SetboUsePigz(boUsePigz=True)
        GwaFile.SetFileHandle(Mode='r')
        GwaFileDCs = GwaFile.ParseToDataContainers()
        GwaFile.Close()
        GwaFile.Cleanup()
        del GwaFile
        LogString = '-- Done ...'
        print LogString
        Log.Write(LogString+'\n')
#       Perform QQ analysis, plot and generate summary
        QQPlotAndSummary(GwaFileDCs,
                         Arguments.QQModes,
                         MtbName,
                         Log)

    if(Arguments.boGeneratePdfReport):
#       Generate summary report in pdf format
        GenerateLaTeXReport(Log,
                            MtbNames)

    #===========================================================================
    # END Do the work!
    #===========================================================================

    #===========================================================================
    # START Finilization
    #===========================================================================
    LogString = '\n**** Done :-)'
    print LogString
    Log.Write(LogString+'\n')
    LogString = Log.GetEndLogString()
    print LogString
    Log.Write(LogString+'\n')
    Log.Close()
    #===========================================================================
    # START Finilization
    #===========================================================================

    return

if(__name__=='__main__'):
    ExecutableName = os.path.abspath(__file__).split('/')[-1]
    main(ExecutableName)