#! /usr/bin/env python
# python modules
import os
import argparse
import re
import scipy
import scipy.stats
import copy
import pylab

# Homebrew modules
import Logger
import ArgumentParser
import File
import DataContainer
import Defines

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
                     'legend.fontsize': 5,
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

def QQPlotAndSummary(DCs=DataContainer.DataContainers,
                     QQModes=argparse.Namespace,
                     Log=Logger):
    boPMode = False
    boSMode = False
    boIMode = False
    if(re.search('P',QQModes)):
        boPMode = True
    if(re.search('S',QQModes)):
        boSMode = True
    if(re.search('I',QQModes)):
        boIMode = True

    # Filter out 'NA' values
    FilterArray  = (DCs.DataContainers['beta'].GetDataArray()!='NA')
    FilterArray *= (DCs.DataContainers['SE'].GetDataArray()!='NA')
    BetaArray    = scipy.real(scipy.compress(FilterArray,DCs.DataContainers['beta'].GetDataArray()).astype(float))
    SEArray      = scipy.real(scipy.compress(FilterArray,DCs.DataContainers['SE'].GetDataArray()).astype(float))

    FilterArray *= (DCs.DataContainers['AF_coded_all'].GetDataArray()!='NA')
    AFArray      = scipy.real(scipy.compress(FilterArray,DCs.DataContainers['AF_coded_all'].GetDataArray()).astype(float))
    CondArray    = [AFArray<=0.5]
    ChoiceArray  = [AFArray]
    MAFArray     = scipy.select(CondArray,ChoiceArray)
    CondArray    = [AFArray>0.5]
    ChoiceArray  = [1.0-AFArray]
    MAFArray    += scipy.select(CondArray,ChoiceArray)
    del CondArray
    del ChoiceArray
    del AFArray
    MAFLevels = copy.deepcopy(Defines.MafLevels)
    MAFLevels.sort()
    ImpLevels = copy.deepcopy(Defines.ImpLevels)
    ImpLevels.sort()

    Chi2Array     = scipy.real(scipy.power((BetaArray/SEArray),2.0)).astype(float)
    FilterArray   = (Chi2Array>=0.0)
    Chi2Array     = scipy.compress(FilterArray,Chi2Array)
    PValArray     = scipy.stats.chi2.cdf(Chi2Array,\
                                         1) # df=1
    LPValArray    = -scipy.log10(PValArray)
    PValExpArray  = scipy.stats.chi2.cdf(scipy.stats.chi2.rvs(1,\
                                                              size=len(PValArray)),\
                                         1)                                          # the 1's are for df=1
    LPValExpArray = -scipy.log10(PValExpArray)

    LogString = '++ Plotting QQ plot in \"p-value\" mode ...'
    print LogString
    Log.Write(LogString+'\n')
    PylabParameters,\
    Rectangle         = PylabGetParams()
    Size = 2.5
    pylab.rcParams.update(PylabParameters)
    PylabFigure = pylab.figure()
    PylabFigure.clf()
    PylabAxis = PylabFigure.add_axes(Rectangle)
    PlotName  = 'Plot.png'
    PylabAxis.scatter(scipy.sort(LPValExpArray),
                      scipy.sort(LPValArray),
                      color=Defines.Colors[0],
                      s=Size)
    for i in range(len(Defines.MafLevels)):
        Color       = Defines.Colors[i+1]
        FilterArray = None
        if(i==0):
            FilterArray = (MAFArray < Defines.MafLevels[i])
        else:
            FilterArray  = (MAFArray >= Defines.MafLevels[i-1])
            FilterArray *= (MAFArray  < Defines.MafLevels[i])
        ObsFP = scipy.sort(scipy.compress(FilterArray,PValArray))
        ObsF  = scipy.sort(scipy.compress(FilterArray,LPValArray))
        scipy.random.shuffle(LPValExpArray)
        ExpF  = scipy.sort(scipy.compress(FilterArray,LPValExpArray))
        PylabAxis.scatter(ExpF,
                          ObsF,
                          color=Color,
                          s=Size)
    MaxLPVal    = LPValArray.max()
    MaxLPValExp = LPValExpArray.max()
    Max         = max(MaxLPVal,MaxLPValExp)+0.5
    PylabAxis.plot([0.0,Max],
                   [0.0,Max],
                   color='black',
                   linestyle='--',
                   linewidth=0.5)
    PylabAxis.set_ylim([0.0,Max])
    PylabAxis.set_xlim([0.0,Max])
    PylabAxis.spines['right'].set_visible(False)
    PylabAxis.spines['top'].set_visible(False)
    PylabAxis.xaxis.set_ticks_position('bottom')
    PylabAxis.yaxis.set_ticks_position('left')
    PylabAxis.grid(True)
    PylabFigure.savefig(PlotName)
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
        QQPlotAndSummary(GwaFileDCs,
                         Arguments.QQModes,
                         Log)
        LogString = '-- Done ...'
        print LogString
        Log.Write(LogString+'\n')

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