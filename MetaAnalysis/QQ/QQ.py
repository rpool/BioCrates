#! /usr/bin/env python
# python modules
import os
import argparse
import re
import scipy
import scipy.stats
import copy

# Homebrew modules
import Logger
import ArgumentParser
import File
import DataContainer
import Defines

def QQPlotAndSummary(DCs=DataContainer.DataContainers,
                     QQModes=argparse.Namespace):
    boPMode = False
    boSMode = False
    if(re.search('P',QQModes)):
        boPMode = True
    if(re.search('S',QQModes)):
        boSMode = True

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

    Chi2Array    = scipy.real(scipy.power((BetaArray/SEArray),2.0)).astype(float)
    FilterArray  = (Chi2Array>=0.0)
    Chi2Array    = scipy.compress(FilterArray,Chi2Array)
    PValArray    = scipy.stats.chi2.cdf(Chi2Array,\
                                        1) # df=1
    LPValArray   = -scipy.log10(PValArray)
    PValExpArray = scipy.stats.chi2.cdf(scipy.stats.chi2.rvs(1,\
                                                             size=len(PValArray)),\
                                        1)                                          # the 1's are for df=1
    LPValExpArray = -scipy.log10(PValExpArray)

#    for i in range(len(LPValArray)):
#        print '@@', LPValExpArray[i], LPValArray[i]
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
        QQPlotAndSummary(GwaFileDCs,
                         Arguments.QQModes)
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