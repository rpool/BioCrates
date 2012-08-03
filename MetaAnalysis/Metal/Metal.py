#! /usr/bin/env python

# python modules
import os
import fnmatch
import sys

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

    MetalHeader  = '#Start metal\n'
    MetalHeader += '#metal\n'
    MetalHeader += '\n'
    MetalHeader  = '#Source this file\n'
    MetalHeader += '#SOURCE thisfile\n'
    MetalHeader += '\n'

    MetalFooter  = '# Analysis\n'
    MetalFooter += 'ANALYZE\n'
    MetalFooter += '\n'
    MetalFooter += '# Exit\n'
    MetalFooter += 'QUIT\n'


    for MtbName in MtbNameFileDCs.DataContainers['0'].GetDataArray():
        ScriptName = MtbName+'_MetalScript.txt'
        fw         = open(ScriptName,'w')
        fw.write(MetalHeader)

        LogString = '++ Generating a \"metal\" script for metabolite\"'+MtbName+'\": \"'+ScriptName+'\" ...'
        print LogString
        Log.Write(LogString+'\n')

        FileNames = fnmatch.filter(Matches,'*_'+MtbName+'_*')

        for CohortName in CohortsFileDCs.DataContainers['0'].GetDataArray():
            FileName = fnmatch.filter(FileNames,
                                      '*'+CohortName+'*')
            if(len(FileName)==1):
                FileName  = FileName[0]
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
            fw.write('#Set Options for describing input files\n')
            fw.write('MARKERLABEL SNPID\n')
            fw.write('PVALUELABEL pval\n')
            fw.write('EFFECTLABEL beta\n')
            fw.write('WEIGHTLABEL n_total\n')
            fw.write('ALLELELABELS coded_all noncoded_all\n')
            fw.write('STRANDLABEL strand_genome\n')
            fw.write('# Strand information\n')
            fw.write('USESTRAND ON\n')
            fw.write('# Column counting per line\n')
            fw.write('COLUMNCOUNTING STRICT\n')
            fw.write('#Process file; assuming that it has been unzipped using e.g. \"gzip -dc\"\n')
            fw.write('PROCESSFILE '+FileName+'\n')

        fw.write('# Outfile\n')
        fw.write('OUTFILE MetaAnalysis_'+MtbName+'_ .tbl\n')
        fw.write('\n')
        fw.write(MetalFooter)
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