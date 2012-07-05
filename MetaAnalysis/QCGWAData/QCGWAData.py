#! /usr/bin/env python
# python modules
import os
import sys
import lxml.etree

# Homebrew modules
import Logger
import ArgumentParser
import File
import Format
import Merge

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

    CommentsPath = os.path.join(os.getcwd(),'Comments')
    if(not os.path.isdir(CommentsPath)):
        os.mkdir(CommentsPath)

    # Perform file format check on extra info file
    LogString = '++ Running file format check on extra info files ...'
    print LogString
    Log.Write(LogString+'\n')
    ExtraInfoFormat = Format.Format()
    ExtraInfoFormat.CheckIfFilesExist(XmlObj=XmlProtocol,
                                      Tag='ExtraInfoFiles',
                                      Log=Log)
    ExtraInfoFormat.SetDelimiter(XmlObj=XmlProtocol,
                                 Log=Log)
    ExtraInfoFormat.SetSplitFunction(Log=Log)
    ExtraInfoFormat.SetColumnFormat(XmlObj=XmlProtocol,
                                    Log=Log)
    ExtraInfoFormat.AppendFilesToExtraInfoFiles(XmlObj=XmlProtocol,
                                                Log=Log)
    boRemoveDuplicateLines = eval(XmlProtocol.getroot().find('Format').find('boRemoveDuplicateLines').text)
    ExtraInfoDCsDict       = ExtraInfoFormat.ParseExtraInfoFiles(Log=Log,
                                                                 boRemoveDuplicateLines=boRemoveDuplicateLines)
    ExtraInfoFormat.CheckFormat(DCsDict=ExtraInfoDCsDict,
                                Log=Log,
                                Path=CommentsPath,
                                FilePreExtName='CheckFormatExtraInfoFiles',
                                FileType='extra info files',
                                XmlObj=XmlProtocol,
                                Tag='ExtraInfoColumns')
    for Key in ExtraInfoDCsDict.iterkeys():
        ExtraInfoDCsDict[Key].DataContainers['SNPID'].SetEntry2IndexDict()
    LogString = '-- Done ...'
    print LogString
    Log.Write(LogString+'\n')

    # Parse MtbNames file and fill MtbNames DataArray
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

        LogString = '++ Processing GWA output files corresponding to \"'+\
                    FileName+\
                    '\" ...'
        print LogString
        Log.Write(LogString+'\n')
        GWADataPath    = XmlProtocol.getroot().find('GWAOutputPath').text.strip()
        GWADataListDir = os.listdir(GWADataPath)
        for N in MtbNameFileDCs.DataContainers['0'].GetDataArray(): # first column, no header
            GWADataFile = None
            Index       = None
            for i in range(len(GWADataListDir)):
                F      = GWADataListDir[i]
                FSplit = F.split('_')
                if(len(FSplit)>3):
                    CompareName = FSplit[2] # if file names are correctly formatted
                    if(CompareName==N):
                        GWADataFile = os.path.join(GWADataPath,F)
                        Index = i
            if(Index!=None):
                del GWADataListDir[Index]

            if(GWADataFile==None):
                LogString  = '  ** Did not find a GWA output file for metabolite \"'+\
                             N+\
                             '\"!\n'
                LogString += '  ** Setting all QC properties to \"False\" or \"Unknown\" for this metabolite!'
                print LogString
                Log.Write(LogString+'\n')
            else:
                LogString  = '  ** Found the GWA output file for metabolite \"'+\
                             N+\
                             '\"!'
                print LogString
                Log.Write(LogString+'\n')
                if(os.path.isfile(GWADataFile) or os.path.islink(GWADataFile)):
                    LogString  = '  ** This file is a regular file or a symbolic link to one!'
                    print LogString
                    Log.Write(LogString+'\n')
                LogString = '  ++ Running file format check on \"'+GWADataFile+'\"  ...'
                print LogString
                Log.Write(LogString+'\n')

                GWAFormat = Format.Format()
                GWAFormat.SetGWADataFileName(Name=GWADataFile)
                GWAFormat.SetDelimiter(XmlObj=XmlProtocol,
                                       Log=Log,
                                       HeadingSpaces='  ')
                GWAFormat.SetSplitFunction(Log=Log,
                                           HeadingSpaces='  ')
                GWAFormat.SetColumnFormat(XmlObj=XmlProtocol,
                                          Log=Log,
                                          HeadingSpaces='  ')
                boRemoveDuplicateLines = eval(XmlProtocol.getroot().find('Format').find('boRemoveDuplicateLines').text)
                GWADCsDict = GWAFormat.ParseGWADataFile(Log=Log,
                                                        HeadingSpaces='  ',
                                                        boRemoveDuplicateLines=boRemoveDuplicateLines)

                GWAFormat.CheckFormat(DCsDict=GWADCsDict,
                                      Log=Log,
                                      HeadingSpaces='  ',
                                      Path=CommentsPath,
                                      FilePreExtName='CheckFormat_'+N,
                                      FileType='GWA data file',
                                      XmlObj=XmlProtocol,
                                      Tag='MtbGWAColumns')
                if(not GWAFormat.GetboColumnFormatOK()):
                    LogString  = '      ** Something is wrong in the column naming!\n'
                    LogString += '      ** Set the proper renaming rules in \"'+Arguments.ProtocolFile+'\"!\n'
                    LogString += '      !! EXITING !!'
                    print LogString
                    Log.Write(LogString+'\n')
                    sys.exit(1)
                boPass = True
                for ColumnId in GWAFormat.GetboFieldFormatOKDict().iterkeys():
                    if(not GWAFormat.GetboFieldFormatOKDict()[ColumnId]):
                        LogString  = '      ** Something is wrong in the field formatting of column \"'+ColumnId+'\"!\n'
                        LogString += '      ** Set the proper renaming rules in \"'+Arguments.ProtocolFile+'\"!'
                        print LogString
                        Log.Write(LogString+'\n')
                        boPass = False
                if(not boPass):
                    LogString = '      !! EXITING !!'
                    print LogString
                    Log.Write(LogString+'\n')
                    sys.exit(1)

                if(len(ExtraInfoFormat.GetExtraInfoFiles())>0):
                    LogString = '    ++ Merging ExtraInfo columns to GWA data file ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    GWADCsDict = Merge.Merge(XmlObj=XmlProtocol,
                                             SourceDCsDict=ExtraInfoDCsDict,
                                             DestDCsDict=GWADCsDict,
                                             SourceColumnTag='SNPID')

                    LogString = '    -- Done ...'
                    print LogString
                    Log.Write(LogString+'\n')

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