#! /usr/bin/env python
# python modules
import os

# Homebrew modules
import Logger
import ArgumentParser
import File
import Format

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
    ExtraInfoFormat.SetDelimiter(XmlObj=XmlProtocol,
                                 Log=Log)
    ExtraInfoFormat.SetSplitFunction(Log=Log)
    ExtraInfoFormat.SetColumnFormat(XmlObj=XmlProtocol,
                                    Log=Log)
    ExtraInfoFormat.AppendFilesToExtraInfoFiles(XmlObj=XmlProtocol,
                                                Log=Log)
    ExtraInfoDCsDict = ExtraInfoFormat.ParseExtraInfoFiles(Log=Log)
    ExtraInfoFormat.CheckFormat(DCsDict=ExtraInfoDCsDict,
                                Log=Log,
                                Path=CommentsPath,
                                FilePreExtName='CheckFormatExtraInfoFiles')
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