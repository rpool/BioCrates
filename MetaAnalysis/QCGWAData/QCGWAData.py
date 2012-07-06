#! /usr/bin/env python
# python modules
import os
import sys
import lxml.etree
import re
import scipy
import scipy.stats

# Homebrew modules
import Logger
import ArgumentParser
import File
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
                    GWAChecksDict[Key]  = Checks.Checks()
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


                for Key in GWADCsDict.iterkeys():
                    MasterFilterArray = None
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

                    LogString = '    ++ QCing all columns: removing rows containing duplicate SNPs (FILTER) ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    if(GWAFiltersDict[Key].GetMaxNDuplicateSNPs()>0):
                        GWADCsDict[Key] = GWAFiltersDict[Key].RemoveDuplicateSNPs(DCs=GWADCsDict[Key])

                        LogString = '      ** Removed '+str(GWAFiltersDict[Key].GetNDeletedDuplicateSNPs())+\
                                    ' row(s) containing duplicate SNPs!'
                        print LogString
                        Log.Write(LogString+'\n')
                    else:
                        LogString = '      ** Removed 0 rows containing duplicate SNPs!'
                        print LogString
                        Log.Write(LogString+'\n')
                    GWAFiltersDict[Key].SetboDuplicateSNPWarning()

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
                    NTotArray  = scipy.copy(GWADCsDict[Key].DataContainers['n_total'].GetDataArray()).astype(int)

#                    PValArray = scipy.stats.t.sf(DataArray,df=NTotArray) does not work!

                    # THIS IS SLOW! MAYBE DO THIS DURING PARSING
                    PValArray = []
                    for i in range(len(DataArray)):
                        T  = DataArray[i]
                        Df = NTotArray[i]
                        PValArray.append(scipy.stats.t.sf(T,\
                                                          Df))
                    PValArray = scipy.array(PValArray)

                    GWADCsDict[Key].DataContainers['PValTTest'].ReplaceDataArray(PValArray.astype(str))

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

                    GWADCsDict[Key].DataContainers['PValWald'].ReplaceDataArray(PValArray.astype(str))
                    del PValArray
                    del DataArray

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

                    LogString = '    -- Done ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    LogString = '    ++ Adding column \"MAF\" to DataContainers ...'
                    print LogString
                    Log.Write(LogString+'\n')

                    GWADCsDict[Key].DataContainers['MAF'] = DataContainer.DataContainer()
                    GWADCsDict[Key].DataContainers['MAF'].SetDataName('MAF')
                    DataArray   = scipy.copy(GWADCsDict[Key].DataContainers['AF_coded_all'].GetDataArray())
                    DataArray   = DataArray.astype(float)
                    FilterArray = DataArray > 0.5
                    DataArray2  = scipy.compress(FilterArray,
                                                 DataArray)
                    DataArray2 -= 0.5
                    scipy.place(DataArray,DataArray>0.5,DataArray2)
                    GWADCsDict[Key].DataContainers['MAF'].ReplaceDataArray(DataArray.astype(str))
                    del DataArray
                    del DataArray2
                    del FilterArray

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