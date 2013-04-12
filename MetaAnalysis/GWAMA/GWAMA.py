#! /usr/bin/env python

# python modules
import os
import fnmatch
import sys
import re

# Homebrew modules
import Logger
import ArgumentParser
import File

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

    InclMtbFileDCs = None
    if(XmlProtocol.getroot().find('InclMtbFile')!=None):
        if(eval(XmlProtocol.getroot().find('InclMtbFile').find('boUse').text)):
            FPath     = XmlProtocol.getroot().find('InclMtbFile').find('Path').text
            FName     = XmlProtocol.getroot().find('InclMtbFile').find('Name').text
            FileName  = os.path.join(FPath,FName)
            LogString = '++ Parsing \"'+FileName+'\" ...'
            print LogString
            Log.Write(LogString+'\n')
            InclMtbFile = File.File(Name=FileName,
                                    boHeader=True)
            InclMtbFile.SetFileHandle(Mode='r')
            boRemoveDuplicateLines = eval(XmlProtocol.getroot().find('Format').find('boRemoveDuplicateLines').text)
            if(boRemoveDuplicateLines):
                NLinesInFile,\
                NLinesInArray  = InclMtbFile.ParseToLineArray()
                LogString = '  ** Removed '+str(NLinesInFile-NLinesInArray)+' duplicate lines!'
                print LogString
                Log.Write(LogString+'\n')
                InclMtbFileDCs = InclMtbFile.LineArray2DataContainers(Delimiter=',')
            else:
                InclMtbFileDCs = InclMtbFile.ParseToDataContainers(Delimiter=',')
            InclMtbFile.Close()
            InclMtbFile.Cleanup()
            del InclMtbFile

            LogString = '-- Done ...'
            print LogString
            Log.Write(LogString+'\n')

    GWAOutputPath = XmlProtocol.getroot().find('GWAOutputPath').text.strip()
    FindREsDict = {}
    for XmlCohort in XmlProtocol.getroot().find('FindRegExps'):
        if(XmlCohort.tag in CohortsFileDCs.DataContainers['0'].GetDataArray()):
            FindREsDict[XmlCohort.tag] = []
        for RE in XmlCohort:
            FindREsDict[XmlCohort.tag].append(RE.text)

    Matches = []
    for Root, DirNames, FileNames in os.walk(GWAOutputPath):
        for Cohort in FindREsDict.iterkeys():
            TmpMatches = []
            for FileName in fnmatch.filter(FileNames,
                                           '*'+Cohort+'*'):
                TmpMatches.append(os.path.join(Root,
                                               FileName))
            for Re in FindREsDict[Cohort]:
                TmpMatches = fnmatch.filter(TmpMatches,
                                            Re)
            Matches.extend(TmpMatches)


    for MtbName in MtbNameFileDCs.DataContainers['0'].GetDataArray():
        ScriptName = MtbName+'_FE_GWAMA.in'
        fw         = open(ScriptName,'w')

        LogString = '++ Generating a \"gwama\" script for metabolite\"'+MtbName+'\": \"'+ScriptName+'\" ...'
        print LogString
        Log.Write(LogString+'\n')

        FileNames = fnmatch.filter(Matches,'*_'+MtbName+'_*')

        for CohortName in CohortsFileDCs.DataContainers['0'].GetDataArray():
            FileName = fnmatch.filter(FileNames,
                                      '*'+CohortName+'*')
            if(len(FileName)==1):
                FileName  = FileName[0]
                Index     = InclMtbFileDCs.DataContainers['Mtb'].GetDataArray().tolist().index(MtbName)
                if(not eval(InclMtbFileDCs.DataContainers[CohortName+'Include[Mtb]'].GetDataArray()[Index])):
                    LogString = '  ** Metabolite \"'+MtbName+'\" ('+CohortName+') is set to be excluded ...'
                    print LogString
                    Log.Write(LogString+'\n')
                    continue
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
            fw.write(re.sub('.gz','',FileName)+'\n')

        fw.close()

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
