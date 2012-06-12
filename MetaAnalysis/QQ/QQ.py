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

def GetPPointsArray(QChi2Array=scipy.array):
    PPointsArray   = [] # probability points array
    ThreeOverEight = 3.0/8.0
    Half           = 0.5
    for Entry in QChi2Array:
        if(Entry <= 10):
            PPointsArray.append((Entry-ThreeOverEight)/(Entry+1-2*ThreeOverEight))
        else:
            PPointsArray.append((Entry-Half)/(Entry+1.0-2.0*Half))
    PPointsArray = scipy.array(PPointsArray)

    return PPointsArray

def LambdaEstimate(PValArray=scipy.array,
                   Filter=True):
    #===========================================================================
    # copied from R-GenABEL.lamdaest
    #===========================================================================
    Estimate   = None
    Ntp        = len(PValArray)

    QChi2Array = scipy.stats.chi2.isf(PValArray,
                                      1)            # quantile function of PValObsArray, df = 1
    if(Filter):
        FilterArray        = (PValArray>=1.0e-8)
        QChi2FilteredArray = scipy.compress(FilterArray,PValArray)
    else:
        QChi2FilteredArray = QChi2Array
    QChi2Array.sort()

    PPointsArray = GetPPointsArray(QChi2FilteredArray)
    PPointsArray = scipy.sort(scipy.stats.chi2.isf((1.0-PPointsArray),1)) # df=1

    FilterArray  = (PPointsArray!=0.0)
    FilterArray *= (QChi2FilteredArray!=0.0)
    FilterArray *= (PPointsArray!=scipy.inf)
    FilterArray *= (QChi2FilteredArray!=scipy.inf)
    FilterArray *= (PPointsArray!=scipy.nan)
    FilterArray *= (QChi2FilteredArray!=scipy.nan)
    FilterArray *= (PPointsArray!=scipy.NaN)
    FilterArray *= (QChi2FilteredArray!=scipy.NaN)
    FilterArray *= (PPointsArray!=scipy.NAN)
    FilterArray *= (QChi2FilteredArray!=scipy.NAN)
    PPointsArray = scipy.compress(FilterArray,PPointsArray)
    QChi2Array   = scipy.compress(FilterArray,QChi2FilteredArray)
#    for i in range(len(PPointsArray)):
#        print '@@',PPointsArray[i],QChi2Array[i]
    P0           = [0.1]
    PBest        = scipy.optimize.leastsq(Residuals,
                                          P0,
                                          args=(QChi2FilteredArray,PPointsArray),
                                          full_output=1,
                                          maxfev=100)
    print PBest
#
##    if(Ntp==1):
##        LambdaEst = 1.0
##    if(PValObsArray.max()<=1.0):
##

    return Estimate

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
    Lambdas    = []
    Lambdas.append(round(LambdaEstimate(PValObsArray),2))
    for i in range(len(MAFLevels)):
        Color       = Defines.Colors[i+1]
        FilterArray = None
        LabelString = None
        if(i==0):
            FilterArray = (MAFArray < MAFLevels[i])
            LabelString = 'MAF < '+str(MAFLevels[i])
            Lambdas.append(0.0)
        else:
            FilterArray  = (MAFArray >= MAFLevels[i-1])
            FilterArray *= (MAFArray  < MAFLevels[i])
            LabelString  = str(MAFLevels[i-1])+r'\leq {\rm MAF} <'+str(MAFLevels[i])
            Lambdas.append(1.0)
        ObsFP = scipy.sort(scipy.compress(FilterArray,PValObsArray))
        ObsF  = scipy.sort(scipy.compress(FilterArray,LPValObsArray))
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

    fw = open(SummaryName,'w')
    print("# Minor Allele Frequency")
    print("# Levels")
    print("Lambdas")
    for i in range(len(MAFLevels)):
        print MAFLevels[i], Lambdas[i]

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
    for i in range(len(ImpQLevels)):
        Color       = Defines.Colors[i+1]
        FilterArray = None
        LabelString = None
        if(i==0):
            FilterArray = (ImpQArray < ImpQLevels[i])
            LabelString = 'ImpQ < '+str(ImpQLevels[i])
        else:
            FilterArray  = (ImpQArray >= ImpQLevels[i-1])
            FilterArray *= (ImpQArray  < ImpQLevels[i])
            LabelString  = str(ImpQLevels[i-1])+r'\leq {\rm ImpQ} <'+str(ImpQLevels[i])
        ObsFP = scipy.sort(scipy.compress(FilterArray,PValObsArray))
        ObsF  = scipy.sort(scipy.compress(FilterArray,LPValObsArray))
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

    return

def PlotQQFilteredOnScore(MtbName=str,
                          PlotPath=str,
                          SummaryPath=str,
                          Chi2ExpArray=scipy.array,
                          Chi2ObsArray=scipy.array,
                          ScoreArray=scipy.array,
                          ScoreLevels=[],
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
    for i in range(len(Defines.MafLevels)):
        Color       = Defines.Colors[i+1]
        FilterArray = None
        LabelString = None
        if(i==0):
            FilterArray = (ScoreArray < Defines.MafLevels[i])
            LabelString = 'Score < '+str(ScoreLevels[i])
        else:
            FilterArray  = (ScoreArray >= ScoreLevels[i-1])
            FilterArray *= (ScoreArray  < ScoreLevels[i])
            LabelString  = str(ScoreLevels[i-1])+r'\leq {\rm Score} <'+str(ScoreLevels[i])
        ObsF  = scipy.sort(scipy.compress(FilterArray,Chi2ObsArray))
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
        AFArray      = scipy.real(scipy.compress(FilterArray,DCs.DataContainers['AF_coded_all'].GetDataArray()).astype(float))
        ImpQArray    = scipy.real(scipy.compress(FilterArray,DCs.DataContainers['oevar_imp'].GetDataArray()).astype(float))
        NTotArray    = scipy.real(scipy.compress(FilterArray,DCs.DataContainers['n_total'].GetDataArray()).astype(float))
        Chi2ObsArray = scipy.real(scipy.power((BetaArray/SEArray),2.0)).astype(float)
        FilterArray  = (Chi2ObsArray>=0.0)
        Chi2ObsArray = scipy.compress(FilterArray,Chi2ObsArray)

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

            PValObsArray  = scipy.stats.chi2.sf(Chi2ObsArray,\
                                                1) # df=1
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
            ScoreArray /= ScoreArray.max()                 # normalize to unity
            ScoreArray *= 100.0                            # normalize to 100%
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
                                  Log)

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
    LogString = '++ Parsing \"'+Arguments.GwaFiles+'\" ...'
    print LogString
    Log.Write(LogString+'\n')
    GwaFiles = File.File(Name=Arguments.GwaFiles,
                         boHeader=False)
    GwaFiles.SetFileHandle(Mode='r')
    GwaFilesDCs = GwaFiles.ParseToDataContainers()
    GwaFiles.Close()
    GwaFiles.Cleanup()
    del GwaFiles
    LogString = '-- Done ...'
    print LogString
    Log.Write(LogString+'\n')

    # Loop over GwaFiles
    for F in GwaFilesDCs.DataContainers['0'].GetDataArray(): # '0' because there should be ONE column && NO header
        LogString = '++ Parsing \"'+F+'\" ...'
        print LogString
        Log.Write(LogString+'\n')
        GwaFile = File.File(Name=F,
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
        MtbName = F.split('_')[2]
        QQPlotAndSummary(GwaFileDCs,
                         Arguments.QQModes,
                         MtbName,
                         Log)

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