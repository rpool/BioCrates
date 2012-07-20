#! /usr/bin/env python
# python modules
import os
import sys
import re
import scipy.stats

# Homebrew modules
import Logger
import ArgumentParser
import File
import HapMap
import Format
import Merge
import Checks
import Filters
import DataContainer

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


    # HapMap stuff
    LogString = '++ Setting up \"HapMap\" data structure ...'
    print LogString
    Log.Write(LogString+'\n')

    HM = HapMap.HapMap()
    HM.ProcessXml(XmlObj=XmlProtocol,
                  Tag='HapMap',
                  Log=Log)
    HM.CheckIfFileExists(Log=Log,
                         HeadingSpaces='')
    HM.ParseDestFile(Log=Log,
                     HeadingSpaces='')
    HM.ProcessLineArray()

    LogString = '-- Done ...'
    print LogString
    Log.Write(LogString+'\n')

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
    Delimiter              = None
    if(XmlProtocol.getroot().find('Format').find('ExtraInfoDelimiter')!=None):
        Delimiter = XmlProtocol.getroot().find('Format').find('ExtraInfoDelimiter').text
    ExtraInfoDCsDict       = ExtraInfoFormat.ParseExtraInfoFiles(Log=Log,
                                                                 boRemoveDuplicateLines=boRemoveDuplicateLines,
                                                                 Delimiter=Delimiter)

    ExtraInfoChecksDict = {}
    for Key in ExtraInfoDCsDict.iterkeys():
        ExtraInfoChecksDict[Key] = Checks.Checks()

    ExtraInfoFormat.CheckFormat(DCsDict=ExtraInfoDCsDict,
                                ChecksDict=ExtraInfoChecksDict,
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

            GWAChecksDict  = {}
            GWAFiltersDict = {}
            if(GWADataFile==None):
                LogString  = '  ** Did not find a GWA output file for metabolite \"'+\
                             N+\
                             '\"!\n'
                LogString += '  ** Setting all QC check outcomes to \"False\" or \"Unknown\" for this metabolite!'
                print LogString
                Log.Write(LogString+'\n')
                GWAChecksDict['GWADataFile']  = Checks.Checks()
                GWAFiltersDict['GWADataFile'] = Filters.Filters()
                GWAChecksDict['GWADataFile'].SetCsvHeader()
                GWAChecksDict['GWADataFile'].WriteCsvHeader('Header.csv')
                GWAChecksDict['GWADataFile'].SetCsvComments('COMMENTS per variable')
                GWAChecksDict['GWADataFile'].WriteCsvComments('Comments.csv')
                GWAChecksDict['GWADataFile'].SetCsvLine(N)
                GWAChecksDict['GWADataFile'].WriteCsvLine('Checks_'+N+'.csv')
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

                for Key in GWADCsDict.iterkeys():
                    GWAChecksDict[Key] = Checks.Checks()
                    GWAChecksDict[Key].SetMtbnameOK(1)
                    GWAChecksDict[Key].SetCsvHeader()
                    GWAChecksDict[Key].WriteCsvHeader('Header.csv')
                    GWAChecksDict[Key].SetCsvComments('COMMENTS per variable')
                    GWAChecksDict[Key].WriteCsvComments('Comments.csv')
                    GWAFiltersDict[Key] = Filters.Filters()

                GWAFormat.CheckFormat(DCsDict=GWADCsDict,
                                      ChecksDict=GWAChecksDict,
                                      Log=Log,
                                      HeadingSpaces='  ',
                                      Path=CommentsPath,
                                      FilePreExtName='CheckFormat_'+N,
                                      FileType='GWA data file',
                                      XmlObj=XmlProtocol,
                                      Tag='MtbGWAColumns')
                GWAFiltersDict[Key].SetMaxNDuplicateSNPs(GWAChecksDict[Key].GetMaxNDuplicateSNPs())

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

                    GWADCsDict = Merge.MergeWithExtraInfo(XmlObj=XmlProtocol,
                                                          SourceDCsDict=ExtraInfoDCsDict,
                                                          DestDCsDict=GWADCsDict,
                                                          SourceColumnTag='SNPID')

                    LogString = '    -- Done ...'
                    print LogString
                    Log.Write(LogString+'\n')

                boDryRun = False
                for XmlTag in XmlProtocol.getroot().find('MtbGWAColumns'):
                    if(XmlTag.find('DryRunFilters')!=None):
                        boDryRun = True
                        break
                if(boDryRun):
                    for Key in GWADCsDict.iterkeys():
                        LogString = '    ++ Dry-running filters ...'
                        print LogString
                        Log.Write(LogString+'\n')

                        LogString = '      ++ QCing column \"n_total\" (DRY-RUN FILTER) ...'
                        print LogString
                        Log.Write(LogString+'\n')

                        GWADCsDict[Key],\
                        FilterTags        = GWAFiltersDict[Key].FilterNTotals(XmlObj=XmlProtocol,
                                                                              DCs=GWADCsDict[Key],
                                                                              ColumnTag='n_total',
                                                                              boDryRun=True)

                        for i in range(len(FilterTags)):
                            FilterTag = FilterTags[i]
                            LogString = '      **'+GWAFiltersDict[Key].GetFilterReportDictDict()['n_total'][FilterTag]
                            LogString = re.sub('\n','\n      **',LogString)
                            if(i<(len(FilterTags)-1)):
                                LogString += '\n'
                            print LogString
                            Log.Write(LogString+'\n')

                        GWAFiltersDict[Key].WriteFilterReport(FileName='DryRunFilter_NTotal_'+N+'.txt',
                                                              Tag='n_total')

                        LogString = '      -- Done ...'
                        print LogString
                        Log.Write(LogString+'\n')

                        LogString = '      ++ QCing column \"SE\" (DRY-RUN FILTER) ...'
                        print LogString
                        Log.Write(LogString+'\n')

                        GWADCsDict[Key],\
                        FilterTags        = GWAFiltersDict[Key].FilterSEs(XmlObj=XmlProtocol,
                                                                          DCs=GWADCsDict[Key],
                                                                          ColumnTag='SE',
                                                                          boDryRun=True)

                        for i in range(len(FilterTags)):
                            FilterTag = FilterTags[i]
                            LogString = '      **'+GWAFiltersDict[Key].GetFilterReportDictDict()['SE'][FilterTag]
                            LogString = re.sub('\n','\n      **',LogString)
                            if(i<(len(FilterTags)-1)):
                                LogString += '\n'
                            print LogString
                            Log.Write(LogString+'\n')

                        GWAFiltersDict[Key].WriteFilterReport(FileName='DryRunFilter_SE_'+N+'.txt',
                                                              Tag='SE')

                        LogString = '      -- Done ...'
                        print LogString
                        Log.Write(LogString+'\n')

                        LogString = '      ++ QCing column \"beta\" (DRY-RUN FILTER) ...'
                        print LogString
                        Log.Write(LogString+'\n')

                        GWADCsDict[Key],\
                        FilterTags        = GWAFiltersDict[Key].FilterBetas(XmlObj=XmlProtocol,
                                                                            DCs=GWADCsDict[Key],
                                                                            ColumnTag='beta',
                                                                            boDryRun=True)

                        for i in range(len(FilterTags)):
                            FilterTag = FilterTags[i]
                            LogString = '      **'+GWAFiltersDict[Key].GetFilterReportDictDict()['beta'][FilterTag]
                            LogString = re.sub('\n','\n      **',LogString)
                            if(i<(len(FilterTags)-1)):
                                LogString += '\n'
                            print LogString
                            Log.Write(LogString+'\n')

                        GWAFiltersDict[Key].WriteFilterReport(FileName='DryRunFilter_Beta_'+N+'.txt',
                                                              Tag='beta')

                        LogString = '      -- Done ...'
                        print LogString
                        Log.Write(LogString+'\n')

                        LogString = '    -- Done ...'
                        print LogString
                        Log.Write(LogString+'\n')


                for Key in GWADCsDict.iterkeys():
                    GWAFiltersDict[Key].SetInitialArrayLength(len(GWADCsDict[Key].DataContainers['SNPID'].GetDataArray()))

                    LogString = '    ++ QCing column \"chr\" (CHECK) ...'
                    print LogString
                    Log.Write(LogString+'\n')
                    if(not GWAChecksDict[Key].CheckChromosomesOK(XmlObj=XmlProtocol,
                                                                 DataArray=GWADCsDict[Key].DataContainers['chr'].GetDataArray(),
                                                                 ColumnTag='chr')):
                        LogString = '      ** Something is wrong in column \"chr\"!'
                        print LogString
                        Log.Write(LogString+'\n')
                    else:
                        LogString = '      ** Column \"chr\" is OK!'
                        print LogString
                        Log.Write(LogString+'\n')
                    LogString = '    -- Done ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    LogString = '    ++ QCing all columns: removing rows containing duplicate SNPs (FILTER) ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    if(GWAFiltersDict[Key].GetMaxNDuplicateSNPs()>0):
                        GWADCsDict[Key] = GWAFiltersDict[Key].RemoveDuplicateSNPs(DCs=GWADCsDict[Key])
                        String = '      ** Removed '+str(GWAFiltersDict[Key].GetNDeletedDuplicateSNPs())+\
                                    ' row(s) containing duplicate SNPs!'
                        GWAFiltersDict[Key].WriteCustomFilterReport(FileName='Filter_SnpIdUnique_'+N+'.txt',
                                                                    String=String)

                        LogString = String
                        print LogString
                        Log.Write(LogString+'\n')
                    else:
                        LogString = '      ** Removed 0 rows containing duplicate SNPs!'
                        print LogString
                        Log.Write(LogString+'\n')
                    GWAFiltersDict[Key].SetboDuplicateSNPWarning()

                    if(GWAFiltersDict[Key].GetboDuplicateSNPWarning()):
                        LogString = '      ** WARNING: REMOVED MORE THAN 100 duplicate SNPs!!'
                        print LogString
                        Log.Write(LogString+'\n')

                    if(len(GWADCsDict[Key].DataContainers['SNPID'].GetDataArray())==0):
                        LogString  = '      !! Data arrays are empty after filtering!\n'
                        LogString += '      !! CONTINUING ...'
                        print LogString
                        Log.Write(LogString+'\n')
                        GWAChecksDict['GWADataFile'].SetCsvLine(N)
                        GWAChecksDict['GWADataFile'].WriteCsvLine('Checks_'+N+'.csv')
                        continue

                    LogString = '    -- Done ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    LogString = '    ++ QCing column \"SE\" (FILTER) ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    GWADCsDict[Key],\
                    FilterTags        = GWAFiltersDict[Key].FilterSEs(XmlObj=XmlProtocol,
                                                                      DCs=GWADCsDict[Key],
                                                                      ColumnTag='SE')

                    for i in range(len(FilterTags)):
                        FilterTag = FilterTags[i]
                        LogString = '      **'+GWAFiltersDict[Key].GetFilterReportDictDict()['SE'][FilterTag]
                        LogString = re.sub('\n','\n      **',LogString)
                        if(i<(len(FilterTags)-1)):
                            LogString += '\n'
                        print LogString
                        Log.Write(LogString+'\n')

                    GWAFiltersDict[Key].WriteFilterReport(FileName='Filter_SE_'+N+'.txt',
                                                          Tag='SE')

                    if(len(GWADCsDict[Key].DataContainers['SNPID'].GetDataArray())==0):
                        LogString  = '      !! Data arrays are empty after filtering!\n'
                        LogString += '      !! CONTINUING ...'
                        print LogString
                        Log.Write(LogString+'\n')
                        GWAChecksDict['GWADataFile'].SetCsvLine(N)
                        GWAChecksDict['GWADataFile'].WriteCsvLine('Checks_'+N+'.csv')
                        continue

                    LogString = '    -- Done ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    LogString = '    ++ QCing column \"beta\" (FILTER) ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    GWADCsDict[Key],\
                    FilterTags        = GWAFiltersDict[Key].FilterBetas(XmlObj=XmlProtocol,
                                                                        DCs=GWADCsDict[Key],
                                                                        ColumnTag='beta')

                    for i in range(len(FilterTags)):
                        FilterTag = FilterTags[i]
                        LogString = '      **'+GWAFiltersDict[Key].GetFilterReportDictDict()['beta'][FilterTag]
                        LogString = re.sub('\n','\n      **',LogString)
                        if(i<(len(FilterTags)-1)):
                            LogString += '\n'
                        print LogString
                        Log.Write(LogString+'\n')

                    GWAFiltersDict[Key].WriteFilterReport(FileName='Filter_Beta_'+N+'.txt',
                                                          Tag='beta')

                    if(len(GWADCsDict[Key].DataContainers['SNPID'].GetDataArray())==0):
                        LogString  = '      !! Data arrays are empty after filtering!\n'
                        LogString += '      !! CONTINUING ...'
                        print LogString
                        Log.Write(LogString+'\n')
                        GWAChecksDict['GWADataFile'].SetCsvLine(N)
                        GWAChecksDict['GWADataFile'].WriteCsvLine('Checks_'+N+'.csv')
                        continue

                    LogString = '    -- Done ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    LogString = '    ++ Mapping GWA data to HapMap build '+\
                                HM.GetDestBuild()+\
                                ' (release '+\
                                HM.GetDestRelease()+\
                                '; '+HM.GetDestRefPanel()+\
                                ') ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    GWADCsDict[Key] = HM.Map(GWADCsDict[Key])

                    RprtName = 'MapHM_'+N+'.txt'
                    LogString = '      ++ Generating report '+RprtName+' in directory \"HapMap\" ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    HM.Report(ReportName=RprtName,
                              SNPIDArray=GWADCsDict[Key].DataContainers['SNPID'].GetDataArray())

                    LogString = '      -- Done ...'
                    print LogString
                    Log.Write(LogString+'\n')
                    HM.Clean()

                    LogString = '    -- Done ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    LogString = '    ++ QCing column \"n_total\" (FILTER) ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    GWADCsDict[Key],\
                    FilterTags        = GWAFiltersDict[Key].FilterNTotals(XmlObj=XmlProtocol,
                                                                          DCs=GWADCsDict[Key],
                                                                          ColumnTag='n_total')

                    for i in range(len(FilterTags)):
                        FilterTag = FilterTags[i]
                        LogString = '      **'+GWAFiltersDict[Key].GetFilterReportDictDict()['n_total'][FilterTag]
                        LogString = re.sub('\n','\n      **',LogString)
                        if(i<(len(FilterTags)-1)):
                            LogString += '\n'
                        print LogString
                        Log.Write(LogString+'\n')

                    GWAFiltersDict[Key].WriteFilterReport(FileName='Filter_NTotal_'+N+'.txt',
                                                          Tag='n_total')

                    if(len(GWADCsDict[Key].DataContainers['SNPID'].GetDataArray())==0):
                        LogString  = '      !! Data arrays are empty after filtering!\n'
                        LogString += '      !! CONTINUING ...'
                        print LogString
                        Log.Write(LogString+'\n')
                        GWAChecksDict['GWADataFile'].SetCsvLine(N)
                        GWAChecksDict['GWADataFile'].WriteCsvLine('Checks_'+N+'.csv')
                        continue

                    LogString = '    -- Done ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    LogString = '    ++ Adding column \"PValTTest\" to DataContainers ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    GWADCsDict[Key].DataContainers['PValTTest'] = DataContainer.DataContainer()
                    GWADCsDict[Key].DataContainers['PValTTest'].SetDataName('PValTtest')
                    DataArray  = scipy.copy(GWADCsDict[Key].DataContainers['beta'].GetDataArray()).astype(float)
                    DataArray /= GWADCsDict[Key].DataContainers['SE'].GetDataArray().astype(float)
                    NTotArray  = scipy.around(GWADCsDict[Key].DataContainers['n_total'].GetDataArray().astype(float)).astype(int)

                    PValArray = map(scipy.stats.t.sf,
                                    scipy.absolute(DataArray).tolist(),
                                    NTotArray.tolist())
                    PValArray = scipy.array(PValArray)
                    GWADCsDict[Key].DataContainers['PValTTest'].ReplaceDataArray(2.0*PValArray)

                    LogString = '    -- Done ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    LogString = '    ++ QCing column \"PValTTest\" against reported column \"pval\" (CHECK) ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    if(not GWAChecksDict[Key].CheckTTestOK(XmlObj=XmlProtocol,
                                                           DataArray=GWADCsDict[Key].DataContainers['pval'].GetDataArray(),
                                                           CheckDataArray=GWADCsDict[Key].DataContainers['PValTTest'].GetDataArray(),
                                                           ColumnTag='pval',
#                                                           boPlot=False,
                                                           boPlot=True,
                                                           MtbName=N)):
                        LogString = '      ** Something is wrong in column \"pval\"!'
                        print LogString
                        Log.Write(LogString+'\n')
                    else:
                        LogString = '      ** Column \"pval\" is OK!'
                        print LogString
                        Log.Write(LogString+'\n')
                    LogString = '    -- Done ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    LogString = '    ++ Adding column \"PValWald\" to DataContainers ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    GWADCsDict[Key].DataContainers['PValWald'] = DataContainer.DataContainer()
                    GWADCsDict[Key].DataContainers['PValWald'].SetDataName('PValWald')
                    DataArray  = scipy.real(scipy.power(DataArray,2.0))
                    PValArray  = scipy.stats.chi2.sf(DataArray,\
                                                     1) # df=1
                    GWADCsDict[Key].DataContainers['PValWald'].ReplaceDataArray(PValArray)
                    del PValArray
                    del DataArray

                    LogString = '    -- Done ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    LogString = '    ++ QCing column \"PValWald\" against reported column \"pval\" (CHECK) ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    if(not GWAChecksDict[Key].CheckNTotalOK(XmlObj=XmlProtocol,
                                                            DataArray=GWADCsDict[Key].DataContainers['pval'].GetDataArray(),
                                                            CheckDataArray=GWADCsDict[Key].DataContainers['PValWald'].GetDataArray(),
                                                            ColumnTag='pval',
#                                                            boPlot=False,
                                                            boPlot=True,
                                                            MtbName=N)):
                        LogString = '      ** Something is wrong in column \"pval\"!'
                        print LogString
                        Log.Write(LogString+'\n')
                    else:
                        LogString = '      ** Column \"pval\" is OK!'
                        print LogString
                        Log.Write(LogString+'\n')
                    LogString = '    -- Done ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    LogString = '    ++ QCing column \"chr\" (FILTER) ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    GWADCsDict[Key],\
                    FilterTags        = GWAFiltersDict[Key].FilterChrs(XmlObj=XmlProtocol,
                                                                       DCs=GWADCsDict[Key],
                                                                       ColumnTag='chr')

                    for i in range(len(FilterTags)):
                        FilterTag = FilterTags[i]
                        LogString = '      **'+GWAFiltersDict[Key].GetFilterReportDictDict()['chr'][FilterTag]
                        LogString = re.sub('\n','\n      **',LogString)
                        if(i<(len(FilterTags)-1)):
                            LogString += '\n'
                        print LogString
                        Log.Write(LogString+'\n')

                    if(len(GWADCsDict[Key].DataContainers['SNPID'].GetDataArray())==0):
                        LogString  = '      !! Data arrays are empty after filtering!\n'
                        LogString += '      !! CONTINUING ...'
                        print LogString
                        Log.Write(LogString+'\n')
                        GWAChecksDict['GWADataFile'].SetCsvLine(N)
                        GWAChecksDict['GWADataFile'].WriteCsvLine('Checks_'+N+'.csv')
                        continue

                    LogString = '    -- Done ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    LogString = '    ++ QCing column \"oevar_imp\" against column \"imputed\" (CHECK) ...'
                    print LogString
                    Log.Write(LogString+'\n')
                    if(not GWAChecksDict[Key].CheckImpQualOK(XmlObj=XmlProtocol,
                                                             ImpQArray=GWADCsDict[Key].DataContainers['oevar_imp'].GetDataArray(),
                                                             ImputedArray=GWADCsDict[Key].DataContainers['imputed'].GetDataArray(),
                                                             ImpQColumnTag='oevar_imp',
                                                             ImputedColumnTag='imputed')):
                        LogString = '      ** Something is wrong in columns \"oevar_imp\" and \"imputed\" !'
                        print LogString
                        Log.Write(LogString+'\n')
                    else:
                        LogString = '      ** Columns \"oevar_imp\" and \"imputed\" are OK!'
                        print LogString
                        Log.Write(LogString+'\n')
                    LogString = '    -- Done ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    LogString = '    ++ QCing column \"n_total\" against column \"imputed\" (CHECK) ...'
                    print LogString
                    Log.Write(LogString+'\n')
                    if(not GWAChecksDict[Key].CheckNTotImpOK(XmlObj=XmlProtocol,
                                                             NTotalArray=GWADCsDict[Key].DataContainers['n_total'].GetDataArray(),
                                                             ImputedArray=GWADCsDict[Key].DataContainers['imputed'].GetDataArray(),
                                                             NTotalColumnTag='n_total',
                                                             ImputedColumnTag='imputed')):
                        LogString = '      ** Something is wrong in columns \"n_total\" and \"imputed\" !'
                        print LogString
                        Log.Write(LogString+'\n')
                    else:
                        LogString = '      ** Columns \"n_total\" and \"imputed\" are OK!'
                        print LogString
                        Log.Write(LogString+'\n')
                    LogString = '    -- Done ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    LogString = '    ++ QCing column \"AF_coded_all\" (FILTER) ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    GWADCsDict[Key],\
                    FilterTags        = GWAFiltersDict[Key].FilterAFCodedAlls(XmlObj=XmlProtocol,
                                                                              DCs=GWADCsDict[Key],
                                                                              ColumnTag='AF_coded_all')

                    for i in range(len(FilterTags)):
                        FilterTag = FilterTags[i]
                        LogString = '      **'+GWAFiltersDict[Key].GetFilterReportDictDict()['AF_coded_all'][FilterTag]
                        LogString = re.sub('\n','\n      **',LogString)
                        if(i<(len(FilterTags)-1)):
                            LogString += '\n'
                        print LogString
                        Log.Write(LogString+'\n')

                    GWAFiltersDict[Key].WriteFilterReport(FileName='Filter_AFCodedAll_'+N+'.txt',
                                                          Tag='AF_coded_all')

                    if(len(GWADCsDict[Key].DataContainers['SNPID'].GetDataArray())==0):
                        LogString  = '      !! Data arrays are empty after filtering!\n'
                        LogString += '      !! CONTINUING ...'
                        print LogString
                        Log.Write(LogString+'\n')
                        GWAChecksDict['GWADataFile'].SetCsvLine(N)
                        GWAChecksDict['GWADataFile'].WriteCsvLine('Checks_'+N+'.csv')
                        continue

                    LogString = '    -- Done ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    LogString = '    ++ Adding column \"MAF\" to DataContainers ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    GWADCsDict[Key].DataContainers['MAF'] = DataContainer.DataContainer()
                    GWADCsDict[Key].DataContainers['MAF'].SetDataName('MAF')
                    DataArray = []
                    for Entry in GWADCsDict[Key].DataContainers['AF_coded_all'].GetDataArray():
                        if(float(Entry)>0.5):
                            DataArray.append(str(1.0-float(Entry)))
                        else:
                            DataArray.append(Entry)
                    DataArray = scipy.array(DataArray)

                    GWADCsDict[Key].DataContainers['MAF'].ReplaceDataArray(DataArray)
                    del DataArray

                    LogString = '    -- Done ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    LogString = '    ++ Adding column \"EMAC\" to DataContainers ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    GWADCsDict[Key].DataContainers['EMAC'] = DataContainer.DataContainer()
                    GWADCsDict[Key].DataContainers['EMAC'].SetDataName('EMAC')
                    DataArray   = scipy.copy(GWADCsDict[Key].DataContainers['n_total'].GetDataArray()).astype(float)
                    DataArray  *= GWADCsDict[Key].DataContainers['MAF'].GetDataArray().astype(float)
                    DataArray  *= GWADCsDict[Key].DataContainers['oevar_imp'].GetDataArray().astype(float)
                    DataArray  *= 2.0
                    GWADCsDict[Key].DataContainers['EMAC'].ReplaceDataArray(DataArray)
                    del DataArray

                    LogString = '    -- Done ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    LogString = '    ++ QCing column \"EMAC\" (FILTER) ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    GWADCsDict[Key],\
                    FilterTags        = GWAFiltersDict[Key].FilterEMACs(XmlObj=XmlProtocol,
                                                                        DCs=GWADCsDict[Key],
                                                                        ColumnTag='EMAC')

                    for i in range(len(FilterTags)):
                        FilterTag = FilterTags[i]
                        LogString = '      **'+GWAFiltersDict[Key].GetFilterReportDictDict()['EMAC'][FilterTag]
                        LogString = re.sub('\n','\n      **',LogString)
                        if(i<(len(FilterTags)-1)):
                            LogString += '\n'
                        print LogString
                        Log.Write(LogString+'\n')

                    GWAFiltersDict[Key].WriteFilterReport(FileName='Filter_EMAC_'+N+'.txt',
                                                          Tag='EMAC')

                    if(len(GWADCsDict[Key].DataContainers['SNPID'].GetDataArray())==0):
                        LogString  = '      !! Data arrays are empty after filtering!\n'
                        LogString += '      !! CONTINUING ...'
                        print LogString
                        Log.Write(LogString+'\n')
                        GWAChecksDict['GWADataFile'].SetCsvLine(N)
                        GWAChecksDict['GWADataFile'].WriteCsvLine('Checks_'+N+'.csv')
                        continue

                    LogString = '    -- Done ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    LogString = '    ++ QCing column \"PValWald\" for inflation (CHECK) ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    MaxNTotal = scipy.copy(GWADCsDict[Key].DataContainers['n_total'].GetDataArray())
                    MaxNTotal = scipy.around(MaxNTotal.astype(float)).astype(int)
                    MaxNTotal = MaxNTotal.max()
                    boOK,\
                    LambdaEst,\
                    SELambdaEst = GWAChecksDict[Key].CheckLambdaOK(XmlObj=XmlProtocol,
                                                                   DataArray=GWADCsDict[Key].DataContainers['PValWald'].GetDataArray(),
                                                                   EMACArray=GWADCsDict[Key].DataContainers['EMAC'].GetDataArray(),
                                                                   MaxNTotal=MaxNTotal,
                                                                   ColumnTag='PValWald',
                                                                   boPlot=True,
                                                                   MtbName=N)
                    if(not boOK):
                        LogString = '      ** Something is wrong in column \"PValWald\"!'
                        print LogString
                        Log.Write(LogString+'\n')
                    else:
                        LogString = '      ** Column \"PValWald\" is OK!'
                        print LogString
                        Log.Write(LogString+'\n')
                    LogString = '      ** Genomic inflation factor Lambda (SE) = '+\
                                str(round(LambdaEst,5))+\
                                ' ('+\
                                str(round(SELambdaEst,5))+\
                                ')'
                    print LogString
                    Log.Write(LogString+'\n')
                    LogString = '    -- Done ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    LogString = '    ++ QCing column \"HapMapMAF\" against column \"MAF\" (CHECK) ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    boOK,\
                    CorrCoeff = GWAChecksDict[Key].CheckScatterFreqsOK(XmlObj=XmlProtocol,
                                                                       HapMapMAFDataArray=GWADCsDict[Key].DataContainers['HapMapMAF'].GetDataArray(),
                                                                       MAFDataArray=GWADCsDict[Key].DataContainers['MAF'].GetDataArray(),
                                                                       ColumnTag='HapMapMAF',
                                                                       boPlot=True,
                                                                       MtbName=N)
                    if(not boOK):
                        LogString = '      ** Something is wrong in columns \"HapMapMAF\" and \"MAF\"!'
                        print LogString
                        Log.Write(LogString+'\n')
                    else:
                        LogString = '      ** Columns \"HapMapMAF\" and \"MAF\" are OK!'
                        print LogString
                        Log.Write(LogString+'\n')
                    LogString = '      **  \"HapMapMAF\" Vs. \"GWAMAF\" correlation coefficient = '+str(round(CorrCoeff,5))
                    print LogString
                    Log.Write(LogString+'\n')

                    GWAChecksDict['GWADataFile'].SetCsvLine(N)
                    GWAChecksDict['GWADataFile'].WriteCsvLine('Checks_'+N+'.csv')

                    LogString = '    -- Done ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    FilteredGWADataPath = os.path.join(GWADataPath,'FILTERED')
                    FilteredGWADataFile = os.path.basename(GWADataFile)
                    if(re.search('.gz',FilteredGWADataFile)):
                        FilteredGWADataFile = re.sub('.gz','',FilteredGWADataFile)
                    HeaderList = ['SNPID',
                                  'chr',
                                  'position',
                                  'coded_all',
                                  'noncoded_all',
                                  'strand_genome',
                                  'beta',
                                  'SE',
                                  'pval',
                                  'AF_coded_all',
                                  'HWE_pval',
                                  'n_total',
                                  'imputed',
                                  'used_for_imp',
                                  'oevar_imp']
                    Header2ColumnDict = {}
                    for Entry in HeaderList:
                        Header2ColumnDict[Entry] = Entry
                    Header2ColumnDict['pval'] = 'PValWald'

                    LogString = '    -- Writing filtered GWA data to \"'+\
                                os.path.join(FilteredGWADataPath,FilteredGWADataFile)+\
                                '\" ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    PValArray  = GWADCsDict[Key].DataContainers['PValWald'].GetDataArray()
                    PValArray  = scipy.around(PValArray,5)
                    TmpArray   = []
                    for Entry in PValArray:
                        TmpArray.append(str(Entry))
                    PValArray = scipy.array(TmpArray)

                    GWADCsDict[Key].DataContainers['PValWald'].ReplaceDataArray(PValArray)

                    GWADCsDict[Key].WriteBioCratesGWAOutput(FileName=FilteredGWADataFile,
                                                            OutPath=FilteredGWADataPath,
                                                            HeaderList=HeaderList,
                                                            Header2ColumnDict=Header2ColumnDict)

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