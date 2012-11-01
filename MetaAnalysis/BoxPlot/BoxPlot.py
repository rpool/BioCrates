#! /usr/bin/env python

# python modules
import os
import fnmatch
import sys
import scipy

# Homebrew modules
import Logger
import ArgumentParser
import File
import Plotting

def main(ExecutableName=str):

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

    # Parse Arguments.ProtocolFile
    LogString = '++ Parsing \"'+Arguments.ProtocolFile+'\" ...'
    print LogString
    Log.Write(LogString+'\n')
    ProtocolFile = File.File(Name=Arguments.ProtocolFile,
                             boHeader=False)
    XmlProtocol  = ProtocolFile.ParseAsXml()
    LogString = '-- Done ...'
    print LogString
    Log.Write(LogString+'\n')

    # Parse NTR scaling file
    if(eval(XmlProtocol.getroot().find('Scaling').find('UseNTRScaling').text)):
        FPath     = XmlProtocol.getroot().find('Scaling').find('Path').text
        FName     = XmlProtocol.getroot().find('Scaling').find('Name').text
        FileName  = os.path.join(FPath,FName)
        LogString = '++ Parsing \"'+FileName+'\" ...'
        print LogString
        Log.Write(LogString+'\n')
        ScalingFile = File.File(Name=FileName,
                                boHeader=True)
        ScalingFile.SetFileHandle(Mode='r')
        boRemoveDuplicateLines = eval(XmlProtocol.getroot().find('Format').find('boRemoveDuplicateLines').text)
        if(boRemoveDuplicateLines):
            NLinesInFile,\
            NLinesInArray  = ScalingFile.ParseToLineArray(Delimiter=',')
            LogString = '  ** Removed '+str(NLinesInFile-NLinesInArray)+' duplicate lines!'
            print LogString
            Log.Write(LogString+'\n')
            ScalingFileDCs = ScalingFile.LineArray2DataContainers()
        else:
            ScalingFileDCs = ScalingFile.ParseToDataContainers()
        ScalingFile.Close()
        ScalingFile.Cleanup()
        del ScalingFile
        LogString = '-- Done ...'
        print LogString
        Log.Write(LogString+'\n')

    # Parse MtbNameFile
    if(eval(XmlProtocol.getroot().find('MtbNameFile').find('boUse').text)):
        FPath     = XmlProtocol.getroot().find('MtbNameFile').find('Path').text
        FName     = XmlProtocol.getroot().find('MtbNameFile').find('Name').text
        FileName  = os.path.join(FPath,FName)
        LogString = '++ Parsing \"'+FileName+'\" ...'
        print LogString
        Log.Write(LogString+'\n')
        MtbNameFile = File.File(Name=FileName,
                                boHeader=False)
        MtbNameFile.SetFileHandle(Mode='r')
        boRemoveDuplicateLines = eval(XmlProtocol.getroot().find('Format').find('boRemoveDuplicateLines').text)
        if(boRemoveDuplicateLines):
            NLinesInFile,\
            NLinesInArray  = MtbNameFile.ParseToLineArray()
            LogString = '  ** Removed '+str(NLinesInFile-NLinesInArray)+' duplicate lines!'
            print LogString
            Log.Write(LogString+'\n')
            MtbNameFileDCs = MtbNameFile.LineArray2DataContainers()
        else:
            MtbNameFileDCs = MtbNameFile.ParseToDataContainers()
        MtbNameFile.Close()
        MtbNameFile.Cleanup()
        del MtbNameFile
        LogString = '-- Done ...'
        print LogString
        Log.Write(LogString+'\n')

    # Parse CohortsFile
    if(eval(XmlProtocol.getroot().find('CohortsFile').find('boUse').text)):
        FPath     = XmlProtocol.getroot().find('CohortsFile').find('Path').text
        FName     = XmlProtocol.getroot().find('CohortsFile').find('Name').text
        FileName  = os.path.join(FPath,FName)
        LogString = '++ Parsing \"'+FileName+'\" ...'
        print LogString
        Log.Write(LogString+'\n')
        CohortsFile = File.File(Name=FileName,
                                boHeader=False)
        CohortsFile.SetFileHandle(Mode='r')
        boRemoveDuplicateLines = eval(XmlProtocol.getroot().find('Format').find('boRemoveDuplicateLines').text)
        if(boRemoveDuplicateLines):
            NLinesInFile,\
            NLinesInArray  = CohortsFile.ParseToLineArray()
            LogString = '  ** Removed '+str(NLinesInFile-NLinesInArray)+' duplicate lines!'
            print LogString
            Log.Write(LogString+'\n')
            CohortsFileDCs = CohortsFile.LineArray2DataContainers()
        else:
            CohortsFileDCs = CohortsFile.ParseToDataContainers()
        CohortsFile.Close()
        CohortsFile.Cleanup()
        del CohortsFile
        LogString = '-- Done ...'
        print LogString
        Log.Write(LogString+'\n')

    GWAOutputPath = XmlProtocol.getroot().find('GWAOutputPath').text.strip()
    FindREs = []
    for RE in XmlProtocol.getroot().find('FindRegExps'):
        FindREs.append(RE.text)

    Matches = []
    for Root, DirNames, FileNames in os.walk(GWAOutputPath):
        for FileName in fnmatch.filter(FileNames,
                                       FindREs[0]):
            Matches.append(os.path.join(Root,
                                        FileName))

    for i in range(1,len(FindREs)):
        RE      = FindREs[i]
        Matches = fnmatch.filter(Matches,
                                 RE)

    for MtbName in MtbNameFileDCs.DataContainers['0'].GetDataArray():
        LogString = '++ Box plotting \"beta\", \"SE\" and \"SD\" as a function of average N for metabolite\"'+MtbName+'\" ...'
        print LogString
        Log.Write(LogString+'\n')

        FileNames = fnmatch.filter(Matches,'*_'+MtbName+'_*')

        CohortGWADCsDict = {}
        CohortList       = []
        for CohortName in CohortsFileDCs.DataContainers['0'].GetDataArray():

            FileName = fnmatch.filter(FileNames,
                                      '*'+CohortName+'*')
            if(len(FileName)==1):
                FileName = FileName[0]
            elif(len(FileName)<1):
                LogString = '  ** No data for metabolite \"'+MtbName+'\" ('+CohortName+') ...'
                print LogString
                Log.Write(LogString+'\n')
                continue
            else:
                LogString  = '  !! Problem in file name matching for metabolite \"'+MtbName+'\" ('+CohortName+') !!\n'
                LogString += '  !! EXITING !!'
                print LogString
                Log.Write(LogString+'\n')
                sys.exit(1)

            CohortList.append(CohortName)

            LogString = '  ++ Parsing \"'+FileName+'\" ...'
            print LogString
            Log.Write(LogString+'\n')

            LogString = '    ++ Constructing data structure for \"'+os.path.basename(FileName)+'\" ...'
            print LogString
            Log.Write(LogString+'\n')

            GWAFile = File.File(Name=FileName,
                                boHeader=True)
            GWAFile.SetboUsePigz(boUsePigz=True)
            GWAFile.SetFileHandle(Mode='r')

            CohortGWADCsDict[CohortName] = GWAFile.ParseTxtFileColumnsToDataContainers(ColumnList=['beta','SE','n_total'])

            GWAFile.Cleanup()

            LogString = '    -- Done ...'
            print LogString
            Log.Write(LogString+'\n')


            LogString = '  -- Done ...'
            print LogString
            Log.Write(LogString+'\n')

#       KORA scaling by NTR summary stats
        KORAStdScalingFactor = None
        if(eval(XmlProtocol.getroot().find('Scaling').find('UseNTRScaling').text)):
            Index                = ScalingFileDCs.DataContainers['Mtb'].GetDataArray().tolist().index(MtbName)
            KORAStdScalingFactor = ScalingFileDCs.DataContainers['Std'].GetDataArray()[Index]
            if(KORAStdScalingFactor!='NA'):
                KORAStdScalingFactor = float(KORAStdScalingFactor)
            else:
                KORAStdScalingFactor = None

#       TwinsUK scaling by 1/log10(e)
        TwinsUKScalingFactor = None
        if(eval(XmlProtocol.getroot().find('Scaling').find('UseTwinsUKScaling').text)):
            TwinsUKScalingFactor = 1.0/scipy.log10(scipy.exp(1.0))

        LogString = '  ++ Plotting \"SE\" box plot ...'
        print LogString
        Log.Write(LogString+'\n')

        XArrayList = []
        YArrayList = []
        NArrayList = []
        for ChrtName in CohortList:
            XArrayList.append(CohortGWADCsDict[ChrtName].DataContainers['n_total'].GetDataArray())
            YArrayList.append(CohortGWADCsDict[ChrtName].DataContainers['SE'].GetDataArray())
            NArrayList.append(CohortGWADCsDict[ChrtName].DataContainers['n_total'].GetDataArray())
        MarkerDict = {}
        MarkerDict['EGCUT']   = 'v'
        MarkerDict['ERF']     = 'o'
        MarkerDict['KORA']    = 's'
        MarkerDict['LLS']     = 'd'
        MarkerDict['NTR']     = 'p'
        MarkerDict['QIMR']    = 'h'
        MarkerDict['TwinsUK'] = '^'

        Plotting.BoxPlotSEPlusConnectingLines(MtbName=MtbName,
                                              DataList=CohortList,
                                              XArrayList=XArrayList,
                                              YArrayList=YArrayList,
                                              MarkerDict=MarkerDict,
                                              KORAScalingFactor=KORAStdScalingFactor,
                                              TwinsUKScaling=TwinsUKScalingFactor)


        LogString = '  -- Done ...'
        print LogString
        Log.Write(LogString+'\n')

        LogString = '  ++ Plotting \"SD\" box plot ...'
        print LogString
        Log.Write(LogString+'\n')

        Plotting.BoxPlotSD(MtbName=MtbName,
                           DataList=CohortList,
                           XArrayList=XArrayList,
                           YArrayList=YArrayList,
                           NArrayList=NArrayList,
                           MarkerDict=MarkerDict,
                           KORAScalingFactor=KORAStdScalingFactor,
                           TwinsUKScaling=TwinsUKScalingFactor)

        LogString = '  -- Done ...'
        print LogString
        Log.Write(LogString+'\n')

        LogString = '  ++ Plotting \"beta\" box plot ...'
        print LogString
        Log.Write(LogString+'\n')

        YArrayList = []
        for ChrtName in CohortList:
            YArrayList.append(CohortGWADCsDict[ChrtName].DataContainers['beta'].GetDataArray())

        Plotting.BoxPlotBeta(MtbName=MtbName,
                             DataList=CohortList,
                             XArrayList=XArrayList,
                             YArrayList=YArrayList,
                             MarkerDict=MarkerDict,
                             KORAScalingFactor=KORAStdScalingFactor,
                             TwinsUKScaling=TwinsUKScalingFactor)

        LogString = '  -- Done ...'
        print LogString
        Log.Write(LogString+'\n')

        LogString = '-- Done ...'
        print LogString
        Log.Write(LogString+'\n')

    #===========================================================================
    # END Do the work!
    #===========================================================================

    #===========================================================================
    # START Finalization
    #===========================================================================
    LogString = '\n**** Done :-)'
    print LogString
    Log.Write(LogString+'\n')
    LogString = Log.GetEndLogString()
    print LogString
    Log.Write(LogString+'\n')
    Log.Close()
    #===========================================================================
    # START Finalization
    #===========================================================================

    return

if(__name__=='__main__'):
    ExecutableName = os.path.abspath(__file__).split('/')[-1]
    main(ExecutableName)