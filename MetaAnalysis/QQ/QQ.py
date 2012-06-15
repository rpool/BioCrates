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

def LinearModel(X=scipy.array,
                Parameters=[]):
    B = Parameters[0]
    Y = B*X
    return Y

def Residuals(Parameters=[],
              X=scipy.array,
              Y=scipy.array):
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
        QChi2Array = scipy.stats.chi2.isf(Array,
                                          1)            # quantile function of PValObsArray, df = 1
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

    P0           = [1.0]
    PBest        = scipy.optimize.leastsq(Residuals,
                                          P0,
                                          args=(PPointsArray,QChi2FilteredArray),
                                          full_output=1,
                                          maxfev=100)
    Estimate     = PBest[0][0]
    # Error estimation of parameter
    Chi2 = scipy.power(PBest[2]['fvec'],2.0).sum()
    Dof  = len(QChi2FilteredArray)-len(P0)-1
    SE   = scipy.real(scipy.sqrt(PBest[1][0,0])*scipy.sqrt(Chi2/float(Dof)))

    return Estimate,SE

def PylabGetParams():
    figwidth_pt   = 246.0 # pt (from revtex \showthe\columnwidth)
    inches_per_pt = 1.0/72.27
    figwidth      = figwidth_pt*inches_per_pt
    golden_mean   = (scipy.sqrt(5.0)-1.0)/2.0 # Aesthetic ratio
    figheight     = figwidth*golden_mean
    fig_size      = [figwidth,figheight]
    params        = {'backend': 'pdf',
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
                     'figure.figsize': fig_size}
    left   = 0.16
    bottom = 0.16
    width  = 0.86-left
    height = 0.95-bottom

    return params,\
           [left, bottom, width, height]

def PlotQQFilteredOnMAF(MtbName=str,
                        PlotPath=str,
                        SummaryPath=str,
                        LPValExpArray=scipy.array,
                        LPValObsArray=scipy.array,
                        PValObsArray=scipy.array,
                        MAFArray=scipy.array,
                        MAFLevels=[],
                        Log=Logger):
    PylabParameters,\
    Rectangle         = PylabGetParams()
    Size              = 2.5
    pylab.rcParams.update(PylabParameters)
    PylabFigure = pylab.figure()
    PylabFigure.clf()
    PylabAxis = PylabFigure.add_axes(Rectangle)
    PlotName  = 'QQPValModeFilteredOnMaf_'+MtbName+'.png'
    PlotName = os.path.join(PlotPath,PlotName)
    PylabAxis.scatter(scipy.sort(LPValExpArray),
                      scipy.sort(LPValObsArray),
                      color=Defines.Colors[0],
                      s=Size,
                      facecolor='None',
                      label=r'${\tt '+MtbName+r'}: {\rm ~all~SNPs}$')
    Lambdas   = []
    SEsLambda = []
    LambdaEst,\
    SELambdaEst = LambdaEstimate(Array=PValObsArray,Filter=False)
    Lambdas.append(LambdaEst)
    SEsLambda.append(SELambdaEst)
    for i in range(len(MAFLevels)):
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
        LambdaEst,\
        SELambdaEst = LambdaEstimate(Array=ObsFP,Filter=False)
        Lambdas.append(LambdaEst)
        SEsLambda.append(SELambdaEst)
        ObsF = scipy.sort(scipy.compress(FilterArray,LPValObsArray))
        scipy.random.shuffle(LPValExpArray)
        ExpF = scipy.sort(scipy.compress(FilterArray,LPValExpArray))
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

    SummaryName = 'QQPValModeFilteredOnMaf_'+MtbName+'.summary.txt'
    SummaryName = os.path.join(SummaryPath,SummaryName)

    LogString = '  ++ Saving MAF filtering summary to \"'+SummaryName+'\" ...'
    print LogString
    Log.Write(LogString+'\n')

    # Set precision
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
                         Log=Logger):
    PylabParameters,\
    Rectangle         = PylabGetParams()
    Size              = 2.5
    pylab.rcParams.update(PylabParameters)
    PylabFigure = pylab.figure()
    PylabFigure.clf()
    PylabAxis = PylabFigure.add_axes(Rectangle)
    PlotName  = 'QQPValModeFilteredOnImpQ_'+MtbName+'.png'
    PlotName = os.path.join(PlotPath,PlotName)
    PylabAxis.scatter(scipy.sort(LPValExpArray),
                      scipy.sort(LPValObsArray),
                      color=Defines.Colors[0],
                      s=Size,
                      facecolor='None',
                      label=r'${\tt '+MtbName+r'}: {\rm ~all~SNPs}$')

    Lambdas   = []
    SEsLambda = []
    LambdaEst,\
    SELambdaEst = LambdaEstimate(Array=PValObsArray,Filter=False)
    Lambdas.append(LambdaEst)
    SEsLambda.append(SELambdaEst)
    for i in range(len(ImpQLevels)):
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
        LambdaEst,\
        SELambdaEst = LambdaEstimate(Array=ObsFP,Filter=False)
        Lambdas.append(LambdaEst)
        SEsLambda.append(SELambdaEst)
        ObsF = scipy.sort(scipy.compress(FilterArray,LPValObsArray))
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


    SummaryName = 'QQPValModeFilteredOnImpQ_'+MtbName+'.summary.txt'
    SummaryName = os.path.join(SummaryPath,SummaryName)
    LogString = '  ++ Saving imputation quality filtering summary to \"'+SummaryName+'\" ...'
    print LogString
    Log.Write(LogString+'\n')

    # Set precision
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
    PylabAxis.scatter(scipy.sort(Chi2ExpArray),
                      scipy.sort(Chi2ObsArray),
                      color=Defines.Colors[0],
                      s=Size,
                      facecolor='None',
                      label=r'${\tt '+MtbName+r'}: {\rm ~all~SNPs}$')

    Lambdas   = []
    SEsLambda = []
    LambdaEst,\
    SELambdaEst = LambdaEstimate(Array=Chi2ObsArray,Filter=True)
    Lambdas.append(LambdaEst)
    SEsLambda.append(SELambdaEst)
    Top20RsIds  = []
    Top20PVals  = []
    Top20RsIds.append([])
    Top20PVals.append([])
    ArgSortArray = scipy.argsort(PValArray)
    for i in range(20):
        Index = ArgSortArray[i]
        Top20RsIds[-1].append(RsIdArray[Index])
        Top20PVals[-1].append(PValArray[Index])
    for i in range(len(ScoreLevels)):
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
        LambdaEst,\
        SELambdaEst = LambdaEstimate(Array=ObsF,Filter=True)
        Lambdas.append(LambdaEst)
        SEsLambda.append(SELambdaEst)
        Top20RsIds.append([])
        Top20PVals.append([])
        RsIdSubArray = scipy.compress(FilterArray,RsIdArray)
        PValSubArray = scipy.compress(FilterArray,PValArray)
        ArgSortArray = scipy.argsort(PValSubArray)
        for i in range(20):
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

    SummaryName = 'QQScoreModeFilteredOnScore_'+MtbName+'.summary.txt'
    SummaryName = os.path.join(SummaryPath,SummaryName)
    LogString = '  ++ Saving score filtering summary to \"'+SummaryName+'\" ...'
    print LogString
    Log.Write(LogString+'\n')

    # Set precision
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
            fw.write(str('# Score level: 0.0 < Score < '+\
                     str(ScoreLevels[i])+'\n'))
        else:
            fw.write(str(ScoreLevels[i-1])+' <= Score < '+\
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
        # set scipy.random.seed
        scipy.random.seed(3168679810)
        # Filter out 'NA' values
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
        PValObsArray  = scipy.stats.chi2.sf(Chi2ObsArray,\
                                            1) # df=1

        # Convert AF to MAF
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
            PValExpArray  = scipy.stats.chi2.sf(scipy.stats.chi2.rvs(1,\
                                                                      size=len(PValObsArray)),\
                                                1,)                                          # the 1's are for df=1
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
                                Log)

            PlotQQFilteredOnImpQ(MtbName,
                                 PlotPath,
                                 SummaryPath,
                                 LPValExpArray,
                                 LPValObsArray,
                                 PValObsArray,
                                 ImpQArray,
                                 ImpQLevels,
                                 Log)

            LogString = '-- Done ...'
            print LogString
            Log.Write(LogString+'\n')
        if(boSMode):
            LogString = '++ Plotting QQ plots in \"score\" mode ...'
            print LogString
            Log.Write(LogString+'\n')

            ScoreArray  = MAFArray*ImpQArray*NTotArray*2.0 # why this expression??
#            ScoreArray /= ScoreArray.max()                 # normalize to unity
#            ScoreArray *= 100.0                            # normalize to 100%
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
    String  = r'\section{QQ $p$--value mode: filtered on MAF}'
    String += '\n\n'
    String += r'Below the QQ plot in $p$--value mode, filtered on '
    String += r'${\rm MAF} \in [0.0,0.5]$, is shown on the left hand side. '
    String += r'The symbol $P$ denotes the $p$--value calculated from the '
    String += r'$\chi^{2}$ distribution of $\beta/SE$ with $df=1$.'
    String += r'For the observed $p$--value $\chi^{2}$ is defined as $\chi^{2}=(\beta/SE)^{2}$. '
    String += r'For the expected $p$--values, a random sample from the $\chi^{2}$ '
    String += r'distribution was taken with the same length as the observed $\chi^{2}$ array.'
    String += r'The right hand side summarizes the estimates for the genomic '
    String += r'inflation factor $\lambda_{\rm est}$ and the associated '
    String += r'standard error ${\rm SE}(\lambda_{\rm est})$, determined at '
    String += r'each MAF level.'
    String += '\n\n'
    String += r'\begin{tabular}{cc}'
    String += '\n'
    File    = re.sub('\.','_',os.path.basename(PlotFile))
    File    = re.sub('_png','.png',File)
    File    = os.path.join('Plots',File)
    os.system('ln -sf '+PlotFile+' '+File)
    String += r'\includegraphics[width=6cm]{../'+File
    String += r'} &'
    String += '\n'
    String += r'\end{tabular}'
    String += '\n\n'
    return String

def LaTeXQQPModeFilteredOnImpQSection(PlotFile=str,
                                      SummaryFile=str):
    String  = r'\section{QQ $p$--value mode: filtered on ImpQ}'
    String += '\n\n'
    String += r'Below the QQ plot in $p$--value mode, filtered on imputation '
    String += r'quality ${\rm ImpQ} \in [0.0,1.0]$, is shown on the left hand side. '
    String += r'The symbol $P$ denotes the $p$--value calculated from the '
    String += r'$\chi^{2}$ distribution $df=1$. '
    String += r'For the observed $p$--value $\chi^{2}$ is defined as $\chi^{2}=(\beta/SE)^{2}$. '
    String += r'For the expected $p$--values, a random sample from the $\chi^{2}$ '
    String += r'distribution was taken with the same length as the observed $\chi^{2}$ array.'
    String += r'The right hand side summarizes the estimates for the genomic '
    String += r'inflation factor $\lambda_{\rm est}$ and the associated '
    String += r'standard error ${\rm SE}(\lambda_{\rm est})$, determined at '
    String += r'each ImpQ level.'
    String += '\n\n'
    String += r'\begin{tabular}{cc}'
    String += '\n'
    File    = re.sub('\.','_',os.path.basename(PlotFile))
    File    = re.sub('_png','.png',File)
    File    = os.path.join('Plots',File)
    os.system('ln -sf '+PlotFile+' '+File)
    String += r'\includegraphics[width=6cm]{../'+File
    String += r'} &'
    String += '\n'
    String += r'\end{tabular}'
    String += '\n\n'
    return String

def LaTeXQQScoreModeFilteredQSection(PlotFile=str,
                                     SummaryFile=str):
    String  = r'\section{QQ $\chi^{2}$ mode: filtered on score}'
    String += '\n\n'
    String += r'Below the QQ plot in $\chi^{2}$ mode, filtered on score '
    String += r'quality ${\rm Score} \in [0.0,100.0]$, is shown on the left hand side. '
    String += r'The observed $\chi^{2}$ values were determined by squaring '
    String += r'$\beta/SE$.'
    String += r'The right hand side summarizes the estimates for the genomic '
    String += r'inflation factor $\lambda_{\rm est}$ and the associated '
    String += r'standard error ${\rm SE}(\lambda_{\rm est})$, determined at '
    String += r'each score level.'
    String += '\n\n'
    String += r'\begin{tabular}{cc}'
    String += '\n'
    File    = re.sub('\.','_',os.path.basename(PlotFile))
    File    = re.sub('_png','.png',File)
    File    = os.path.join('Plots',File)
    os.system('ln -sf '+PlotFile+' '+File)
    String += r'\includegraphics[width=6cm]{../'+File
    String += r'} &'
    String += '\n'
    String += r'\end{tabular}'
    String += '\n\n'
    return String

def LaTeXPreamble():
    String  =r'\documentclass[pre,amsmath,onecolumn,floatfix,fleqn,a4paper,superscriptaddress]{revtex4}'
    String += '\n'
    String +=r'\usepackage{latexsym}'
    String += '\n'
    String +=r'\usepackage{amsmath,amsfonts,amssymb}'
    String += '\n'
    String +=r'\usepackage[dvips]{graphicx}'
    String += '\n'
    String +=r'\usepackage{subfigure}'
    String += '\n'
    String +=r'\usepackage{float}'
    String += '\n'
    String +=r'\usepackage{pifont}'
    String += '\n'
    String +=r'\usepackage{color}'
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
    String  = r'\title{QQ Analysis Report for Metabolite '
    String += r'{\tt '
    String += MtbName
    String += r'}'
    String += r'}'
    String += '\n'
    return String

def LaTeXDate():
    String  = r'\date{\today}'
    String += '\n'
    String += '\n'
    return String

def LaTeXBeginDocument():
    String  = r'\begin{document}'
    String += '\n'
    return String

def LaTeXMakeTitle():
    String  = r'\maketitle'
    String += '\n'
    String += '\n'
    return String

def LaTeXEndDocument():
    String  = r'\end{document}'
    String += '\n'
    return String


def GenerateLaTeXReport(Log=Logger,
                        MtbNames=[]):
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
        del GwaDataFileList[DelIndex]
        LogString = '++ Parsing \"'+GwaFileName+'\" ...'
        print LogString
        Log.Write(LogString+'\n')
#        GwaFile = File.File(Name=GwaFileName,
#                            boHeader=True)
#        GwaFile.SetboUsePigz(boUsePigz=True)
#        GwaFile.SetFileHandle(Mode='r')
#        GwaFileDCs = GwaFile.ParseToDataContainers()
#        GwaFile.Close()
#        GwaFile.Cleanup()
#        del GwaFile
        LogString = '-- Done ...'
        print LogString
        Log.Write(LogString+'\n')
#        QQPlotAndSummary(GwaFileDCs,
#                         Arguments.QQModes,
#                         MtbName,
#                         Log)

    if(Arguments.boGeneratePdfReport):
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