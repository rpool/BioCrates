#! /usr/bin/env python
import argparse
import os
import sys
import re
import gzip
import itertools
import scipy
import scipy.stats

import SummarizeFigsUsingTeX
import BioCratesAnalyticalRanges
import Plotting
import ExcelParser
import ArgumentParser
import QC
import Logger

def main(ExecutableName):
    ArgParser,\
    Arguments   = ArgumentParser.ParseArguments()

    SampleDataDict = None
    if(not Arguments.boQCDone):
        Ext = 'qc'
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
#       Do QC and Impute
        SampleDataDict = QCAndImpute(Arguments,
                                     Log)
        LogString = '**** Done :-)'
        print LogString
        Log.Write(LogString+'\n')
        LogString = Log.GetEndLogString()
        print LogString
        Log.Write(LogString+'\n')
        Log.Close()

    Ext = 'postprocess'
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
#   Perform postprocessing procedure
    PostProcess(SampleDataDict,
                Arguments,
                Log)
    LogString = '**** Done :-)'
    print LogString
    Log.Write(LogString+'\n')
    LogString = Log.GetEndLogString()
    print LogString
    Log.Write(LogString+'\n')
    Log.Close()

    return

def PostProcess(SampleDataDict=None,
                Arguments=argparse.Namespace,
                Log=Logger):
    if(SampleDataDict==None):
        LogString = '**** Parsing Excel book \"'+Arguments.ExcelSampleDataFileName+'\" containing the sample data ...'
        print LogString
        Log.Write(LogString+'\n')
        SampleDataDict = ExcelParser.ParseExcellBook(Arguments.ExcelSampleDataFileName,
                                                     Arguments.ExcelSampleSheetName,
                                                     Arguments.boRemoveCustomMetabolites,
                                                     Log)
    FileList = []
    FileList.append('MiceImpute.report.gz')
    FileList.append('QCBioCratesMetabolomics.py.qc.log')
#   check if timestamps are the same
    TimeStamps = []
    for File in FileList:
        fr = None
        if(re.search('.gz',File)):
            fr = gzip.open(File,'rb')
        else:
            fr = open(File,'r')
        boFoundTag = False
        TimeStamp  = ''
        for Line in fr.readlines():
            if(Line.strip()=='## END TIMESTAMP'):
                break
            if(boFoundTag):
                TimeStamp += Line
            if(Line.strip()=='## START TIMESTAMP'):
                boFoundTag = True
        fr.close()
        TimeStamps.append(TimeStamp)
    TimeStamps = list(set(TimeStamps))
    if(len(TimeStamps)>1):
        LogString = '** TimeStamps of input files differ ...\n'
        LogString = '** EXITING ...'
        print LogString
        Log.Write(LogString+'\n')
        sys.exit(1)

    fr   = gzip.open('MiceImpute.report.gz','rb')
    FMem = fr.readlines()
    fr.close()

    boFoundCompleteTag = False
    boFoundDataTag     = False
    Data               = ''
    Complete           = ''
    for Line in FMem:
        if(Line.strip()=='## START DATA'):
            boFoundDataTag = True
            continue
        if(Line.strip()=='## END DATA'):
            boFoundDataTag = False
            continue
        if(Line.strip()=='## START COMPLETE'):
            boFoundCompleteTag = True
            continue
        if(Line.strip()=='## END COMPLETE'):
            boFoundCompleteTag = False
            continue
        if(boFoundCompleteTag):
            Complete += Line
        if(boFoundDataTag):
            Data += Line

    Data       = Data.split('\n')
    DataDict   = {}
    DataHeader = Data[0].strip().split()
    DataHeader.reverse()
    for i in range(len(DataHeader)):
        Entry = int(re.sub('X','',DataHeader[i]))
        if(not DataDict.has_key(Entry)):
            DataDict[Entry] = []
        for j in range(1,len(Data)):
            Line = Data[j].strip().split()
            if(len(Line)==0):
                continue
            Value = Line[-1-i]
            if(Value!='NA'):
                Value = float(Value)
            DataDict[Entry].append(Value)

    Complete       = Complete.split('\n')
    CompleteDict   = {}
    CompleteHeader = Complete[0].strip().split()
    CompleteHeader.reverse()
    for i in range(len(CompleteHeader)):
        Entry = int(re.sub('X','',CompleteHeader[i]))
        if(not CompleteDict.has_key(Entry)):
            CompleteDict[Entry] = []
        for j in range(1,len(Complete)):
            Line = Complete[j].strip().split()
            if(len(Line)==0):
                continue
            Value = Line[-1-i]
            if(Value!='NA'):
                Value = float(Value)
            CompleteDict[Entry].append(Value)
    del FMem

    for Key,Value in DataDict.iteritems():
        SampleDataDict[Key].SetQCedDataArray(Value)
    for Key,Value in CompleteDict.iteritems():
        SampleDataDict[Key].SetImputedDataArray(Value)

    if(Arguments.boPlotDistributions):
        LogString = '**** Parsing analytical ranges from \"'+str(Arguments.AnalyticalRangesExcelFileName)+'\" ...'
        print LogString
        Log.Write(LogString+'\n')
        AnalRanges = BioCratesAnalyticalRanges.BioCratesAnalyticalRanges(Arguments.AnalyticalRangesExcelFileName)
        LogString = '**** Plotting sample data distributions ...'
        print LogString
        Log.Write(LogString+'\n')
#        PlotNames  = Plotting.PlotDistributions(SampleDataDict,
#                                                AnalRanges,
#                                                'SAMPLE'+Arguments.PlotBaseName,
#                                                '.png',
#                                                300,
#                                                Log)
#        SummarizeFigsUsingTeX.SummarizePlots(PlotNames,
#                                             'SAMPLE'+Arguments.SummaryPlotBaseName.split('.')[0],
#                                             '.png',
#                                             Log)
        LnSpacePlotNames  = Plotting.PlotLnSpaceDistributions(SampleDataDict,
                                                              AnalRanges,
                                                              'SAMPLE_LNSPACE'+Arguments.PlotBaseName,
                                                              '.png',
                                                              300,
                                                              Log)
        SummarizeFigsUsingTeX.SummarizePlots(LnSpacePlotNames,
                                             'SAMPLE_LNSPACE'+Arguments.SummaryPlotBaseName.split('.')[0],
                                             '.png',
                                             Log)

    if(Arguments.boCorrelateCombinations):
        ImputedMetabolites = []
        for Key in SampleDataDict.iterkeys():
            if(SampleDataDict[Key].GetImputedDataArray()!=None):
                ImputedMetabolites.append(Key)
        ImputedCombinations = list(itertools.combinations(range(len(ImputedMetabolites)),2))

        fw = open('CombinationCorrelations.txt','w')
        fw.write('#')
        FormatString = '{0:>7s}'
        fw.write(FormatString.format('IndexA'))
        fw.write(FormatString.format('IndexB'))
        FormatString = '{0:>16s}'
        fw.write(FormatString.format('NameA'))
        fw.write(FormatString.format('NameB'))
        FormatString = '{0:>14s}'
        fw.write(FormatString.format('a'))
        fw.write(FormatString.format('c'))
        fw.write(FormatString.format('R^2'))
        fw.write(FormatString.format('p-value'))
        fw.write(FormatString.format('stderr'))
        fw.write(' (of regression ln[B] = a ln[A] + c)')
        fw.write('\n')
        for Combination in ImputedCombinations:
            CombinationDataArray = []
            CombinationNames     = []
            CombinationIndices   = []
            for MetaboliteKey in Combination:
                CombinationDataArray.append(SampleDataDict[ImputedMetabolites[MetaboliteKey]].GetImputedDataArray())
                CombinationNames.append(SampleDataDict[ImputedMetabolites[MetaboliteKey]].GetMetaboliteConventionName())
                CombinationIndices.append(SampleDataDict[ImputedMetabolites[MetaboliteKey]].GetMetaboliteIndex())
#            Slope,\
#            Intercept,\
#            RValue,\
#            PValue,\
#            StdErr       = scipy.stats.linregress(scipy.array(CombinationDataArray[0]),
#                                                  scipy.array(CombinationDataArray[1]))
            Slope,\
            Intercept,\
            RValue,\
            PValue,\
            StdErr       = scipy.stats.linregress(scipy.log(scipy.array(CombinationDataArray[0])+1.0),
                                                  scipy.log(scipy.array(CombinationDataArray[1])+1.0))
            FormatString = '{0:>8s}'
            for Entry in CombinationIndices:
                fw.write(FormatString.format(str(Entry)))
                FormatString = '{0:>7s}'
            FormatString = '{0:>16s}'
            for Entry in CombinationNames:
                fw.write(FormatString.format(str(Entry)))
            FormatString = '{0:>+14.5e}'
            fw.write(FormatString.format(Slope))
            fw.write(FormatString.format(Intercept))
            fw.write(FormatString.format(float(scipy.real(scipy.power(RValue,2.0)))))
            fw.write(FormatString.format(PValue))
            fw.write(FormatString.format(StdErr))
            fw.write('\n')

#            Slope,\
#            Intercept,\
#            RValue,\
#            PValue,\
#            StdErr       = scipy.stats.linregress(scipy.array(CombinationDataArray[1]),
#                                                  scipy.array(CombinationDataArray[0]))
            Slope,\
            Intercept,\
            RValue,\
            PValue,\
            StdErr       = scipy.stats.linregress(scipy.log(scipy.array(CombinationDataArray[1])+1.0),
                                                  scipy.log(scipy.array(CombinationDataArray[0])+1.0))
            FormatString = '{0:>8s}'
            CombinationIndices.reverse()
            for Entry in CombinationIndices:
                fw.write(FormatString.format(str(Entry)))
                FormatString = '{0:>7s}'
            FormatString = '{0:>16s}'
            CombinationNames.reverse()
            for Entry in CombinationNames:
                fw.write(FormatString.format(str(Entry)))
            FormatString = '{0:>+14.5e}'
            fw.write(FormatString.format(Slope))
            fw.write(FormatString.format(Intercept))
            fw.write(FormatString.format(float(scipy.real(scipy.power(RValue,2.0)))))
            fw.write(FormatString.format(PValue))
            fw.write(FormatString.format(StdErr))
            fw.write('\n')

        fw.close()

    return

def QCAndImpute(Arguments=argparse.Namespace,
                Log=Logger):

    LogString = '**** Parsing Excel book \"'+Arguments.ExcelSampleDataFileName+'\" containing the sample data ...'
    print LogString
    Log.Write(LogString+'\n')
    SampleDataDict = ExcelParser.ParseExcellBook(Arguments.ExcelSampleDataFileName,
                                                 Arguments.ExcelSampleSheetName,
                                                 Arguments.boRemoveCustomMetabolites,
                                                 Log)
    LogString = '**** Parsing Excel book \"'+Arguments.ExcelReferenceDataFileName+'\" containing the reference data ...'
    print LogString
    Log.Write(LogString+'\n')
    ReferenceDataDict = ExcelParser.ParseExcellBook(Arguments.ExcelReferenceDataFileName,
                                                    Arguments.ExcelReferenceSheetName,
                                                    Arguments.boRemoveCustomMetabolites,
                                                    Log)
    LogString = '**** Parsing analytical ranges from \"'+str(Arguments.AnalyticalRangesExcelFileName)+'\" ...'
    print LogString
    Log.Write(LogString+'\n')
    AnalRanges = BioCratesAnalyticalRanges.BioCratesAnalyticalRanges(Arguments.AnalyticalRangesExcelFileName)
    for Value in SampleDataDict.itervalues():
        Value.SetLODFromDocumentation(AnalRanges)
        Value.SetLLOQFromDocumentation(AnalRanges)
        Value.SetULOQFromDocumentation(AnalRanges)
    for Value in ReferenceDataDict.itervalues():
        Value.SetLODFromDocumentation(AnalRanges)
        Value.SetLLOQFromDocumentation(AnalRanges)
        Value.SetULOQFromDocumentation(AnalRanges)

    if(Arguments.boPlotDistributions):
        LogString = '**** Plotting sample data distributions ...'
        print LogString
        Log.Write(LogString+'\n')
        PlotNames  = Plotting.PlotDistributions(SampleDataDict,
                                                AnalRanges,
                                                'SAMPLE'+Arguments.PlotBaseName,
                                                '.png',
                                                300,
                                                Log)
        SummarizeFigsUsingTeX.SummarizePlots(PlotNames,
                                             'SAMPLE'+Arguments.SummaryPlotBaseName.split('.')[0],
                                             '.png',
                                             Log)
        LogString = '**** Plotting reference data distributions ...'
        print LogString
        Log.Write(LogString+'\n')
        PlotNames  = Plotting.PlotDistributions(ReferenceDataDict,
                                                AnalRanges,
                                                'REFERENCE'+Arguments.PlotBaseName,
                                                '.png',
                                                300,
                                                Log)
        SummarizeFigsUsingTeX.SummarizePlots(PlotNames,
                                             'REFERENCE'+Arguments.SummaryPlotBaseName.split('.')[0],
                                             '.png',
                                             Log)

    ReferenceQCMetaboliteContainers = []
    SampleQCMetaboliteContainers    = []
    if(Arguments.boQCReference):
        LogString = '**** Performing quality control using the reference data ...'
        print LogString
        Log.Write(LogString+'\n')
        ReferenceQCMetaboliteContainers = QC.QCReference(ReferenceDataDict,
                                                         Log)
    if(Arguments.boQCSample):
        LogString = '**** Performing quality control using the sample data ...'
        print LogString
        Log.Write(LogString+'\n')
        SampleQCMetaboliteContainers = QC.QCSample(SampleDataDict,
                                                   ReferenceQCMetaboliteContainers,
                                                   Arguments.DLHandling,
                                                   Log)
    if(Arguments.boQCReference or Arguments.boQCSample):
        if(Arguments.boImpute):
            LogString = '**** Processing sample and/or reference QC results ...'
            print LogString
            Log.Write(LogString+'\n')
            fw = open('QC.report','w')
            fw.write('## START TIMESTAMP\n')
            fw.write(str(Log.GetStartDate())+'\n')
            fw.write('## END TIMESTAMP\n')
            MetaboliteNameExclusionList = []
            SampleIDExclusionList       = []
            ProcessedDataDict,\
            NameExclusionList,\
            IDExclusionList     = QC.ProcessDataDict(ReferenceQCMetaboliteContainers,
                                                     SampleQCMetaboliteContainers,
                                                     SampleDataDict)
            MetaboliteNameExclusionList.extend(NameExclusionList)
            SampleIDExclusionList.extend(IDExclusionList)
            ProcessedDataDict,\
            NameExclusionList    = QC.QCFinal(ProcessedDataDict,
                                              Log)
            MetaboliteNameExclusionList.extend(NameExclusionList)
            MetaboliteNameExclusionList = list(set(MetaboliteNameExclusionList))
            MetaboliteNameExclusionList.sort()
            fw.write('## START EXCLUDEDMETABOLITES\n')
            for i in range(len(MetaboliteNameExclusionList)):
                Entry = MetaboliteNameExclusionList[i]
                fw.write('{0:>4}'.format(str(i))+'{0:>25}'.format(str(Entry)))
                Entry2 = None
                for Value in SampleDataDict.itervalues():
                    if(Value.GetMetaboliteName()==Entry):
                        Entry2 = Value.GetMetaboliteConventionName()
                fw.write('{0:>25}'.format(Entry2))
                fw.write('\n')
            fw.write('## END EXCLUDEDMETABOLITES\n')
            fw.write('## START EXCLUDEDSAMPLEIDS\n')
            SampleIDExclusionList.sort()
            for i in range(len(SampleIDExclusionList)):
                Entry = SampleIDExclusionList[i]
                fw.write(str(i)+' '+str(Entry)+'\n')
            fw.write('## END EXCLUDEDSAMPLEIDS\n')
            fw.close()
            Impute = QC.RMiceImpute()
            Impute.Impute(Arguments,
                          ProcessedDataDict,
                          Log)

#   Now we have to combine the datacontainers and perform a last QC run (Missing Value rate>5%)
#   Then impute the missings using mice (R)
    if(not(Arguments.boPlotDistributions) and
       not(Arguments.boQCReference) and
       not(Arguments.boQCSample)):
        LogString = '**** Nothing done ...'
        print LogString
        Log.Write(LogString+'\n')

    return SampleDataDict

if(__name__=='__main__'):

    ExecutableName = os.path.abspath(__file__).split('/')[-1]
    main(ExecutableName)
