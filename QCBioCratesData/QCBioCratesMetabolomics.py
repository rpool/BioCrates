#! /usr/bin/env python
import argparse
import os
import sys
import re
import gzip
import itertools
import scipy
import scipy.stats
import copy

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
#   Remove excluded samples
    fr       = open('QC.report','r')
    QCReport = fr.readlines()
    fr.close()
    ExcludedSamples = []
    for i in range(QCReport.index('## START EXCLUDEDSAMPLEIDS\n')+1,QCReport.index('## END EXCLUDEDSAMPLEIDS\n')):
        ExcludedSamples.append(QCReport[i].strip().split()[-1])

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
    DataNames = []
    for i in range(len(DataHeader)):
        Entry = int(re.sub('X','',DataHeader[i]))
        DataNames.append(SampleDataDict[Entry].GetMetaboliteConventionName())
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

    CompleteDataNames = []
    for Key in CompleteDict.iterkeys():
        CompleteDataNames.append(SampleDataDict[Key].GetMetaboliteConventionName())

    for Key,Value in DataDict.iteritems():
        SampleDataDict[Key].SetQCedDataArray(Value)
    for Key,Value in CompleteDict.iteritems():
        SampleDataDict[Key].SetImputedDataArray(Value)

    SIDKey = None
    for Key in SampleDataDict.iterkeys():
        if(SampleDataDict[Key].GetDataName()=='Sample Identification'):
            SIDKey = Key
            break
    SampleIDList = copy.deepcopy(SampleDataDict[SIDKey].GetDataArray())
    DelIndexList = []
    for S in ExcludedSamples:
        DelIndexList.append(SampleIDList.index(S))
    DelIndexList.sort()
    DelIndexList.reverse()
    for i in DelIndexList:
        del SampleIDList[i]

    LogString = '**** Writing QCed data to \"QCedData.csv\" ...'
    print LogString
    Log.Write(LogString+'\n')
    fw = open('QCedData.csv','w')
    String = ['']
    String.append('sample.id')
    for Key in sorted(DataDict.keys()):
        String.append(SampleDataDict[Key].GetMetaboliteConventionName())
    fw.write(','.join(String)+'\n')
    for i in range(len(DataDict[DataDict.keys()[0]])):
        String = []
        String.append(str(i+1))
        String.append(SampleIDList[i])
        for j in sorted(DataDict.keys()):
            String.append(str(DataDict[j][i]))
        fw.write(','.join(String)+'\n')
    fw.close()

    LogString = '**** Writing ln-transformed QCed data to \"LnQCedData.csv\" ...'
    print LogString
    Log.Write(LogString+'\n')
    fw = open('LnQCedData.csv','w')
    String = ['']
    String.append('sample.id')
    for Key in sorted(DataDict.keys()):
        String.append(SampleDataDict[Key].GetMetaboliteConventionName())
    fw.write(','.join(String)+'\n')
    for i in range(len(DataDict[DataDict.keys()[0]])):
        String = []
        String.append(str(i+1))
        String.append(SampleIDList[i])
        for j in sorted(DataDict.keys()):
            if(type(DataDict[j][i])==str):
                String.append(DataDict[j][i])
            elif(DataDict[j][i]<=0.0):
                String.append('NA')
            else:
                String.append(str(scipy.log(DataDict[j][i])))
        fw.write(','.join(String)+'\n')
    fw.close()

    LogString = '**** Writing QCed data to \"QCedAndImputedData.csv\" ...'
    print LogString
    Log.Write(LogString+'\n')
    fw = open('QCedAndImputedData.csv','w')
    String = ['']
    String.append('sample.id')
    for Key in sorted(CompleteDict.keys()):
        String.append(SampleDataDict[Key].GetMetaboliteConventionName())
    fw.write(','.join(String)+'\n')
    for i in range(len(CompleteDict[CompleteDict.keys()[0]])):
        String = []
        String.append(str(i+1))
        String.append(SampleIDList[i])
        for j in sorted(CompleteDict.keys()):
            String.append(str(CompleteDict[j][i]))
        fw.write(','.join(String)+'\n')
    fw.close()

    LogString = '**** Writing ln-transformed QCed data to \"LnQCedAndImputedData.csv\" ...'
    print LogString
    Log.Write(LogString+'\n')
    fw = open('LnQCedAndImputedData.csv','w')
    String = ['']
    String.append('sample.id')
    for Key in sorted(CompleteDict.keys()):
        String.append(SampleDataDict[Key].GetMetaboliteConventionName())
    fw.write(','.join(String)+'\n')
    for i in range(len(CompleteDict[CompleteDict.keys()[0]])):
        String = []
        String.append(str(i+1))
        String.append(SampleIDList[i])
        for j in sorted(CompleteDict.keys()):
            if(type(CompleteDict[j][i])==str):
                String.append(CompleteDict[j][i])
            elif(CompleteDict[j][i]<=0.0):
                String.append('NA')
            else:
                String.append(str(scipy.log(CompleteDict[j][i])))
        fw.write(','.join(String)+'\n')
    fw.close()

    if(Arguments.boNormalityTest):
        ExcludedMtbList = []
        fr              = open('QC.report','r')
        boStart         = False
        for Line in fr:
            if(Line.strip()=='## END EXCLUDEDMETABOLITES'):
                break
            if(Line.strip()=='## START EXCLUDEDMETABOLITES'):
                boStart = True
                continue
            if(boStart):
                Name = Line.strip().split()
                Name = ' '.join(Name[1:-1])
                ExcludedMtbList.append(Name)
        fr.close()
        NormalityTest(SampleDataDict,
                      ExcludedMtbList,
                      Log)

    if(Arguments.boPlotDistributions):
        LogString = '**** Parsing analytical ranges from \"'+str(Arguments.AnalyticalRangesExcelFileName)+'\" ...'
        print LogString
        Log.Write(LogString+'\n')
        AnalRanges = BioCratesAnalyticalRanges.BioCratesAnalyticalRanges(Arguments.AnalyticalRangesExcelFileName)
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

def NormalityTest(DataDict=dict,
                  ExcludedMtbList=list,
                  Log=Logger):
    LogString  = '**** Performing normality tests on the concentrations, '
    LogString += 'the ln-transformed concentrations and all possible ratios in ln-space'
    print LogString
    Log.Write(LogString+'\n')

    LogString = '  ** Performing normality tests on the concentrations ...'
    print LogString
    Log.Write(LogString+'\n')

    fw = open('NormalityCheckConcentrations.csv','w')
    fw.write('metabolite,KS test statistic, KS p-value\n')
    for Value in DataDict.itervalues():
        if(Value.GetMetaboliteName() and
           (not Value.GetMetaboliteName() in ExcludedMtbList)):
            fw.write(Value.GetMetaboliteConventionName()+',')

            DataArray = None
            if(Value.GetImputedDataArray()!=None):
                DataArray = scipy.real(scipy.array(Value.GetImputedDataArray()))
            else:
                DataArray  = scipy.real(scipy.array(Value.GetDataArray()))
            Mu         = DataArray.mean()
            Sd         = DataArray.std(ddof=1)
            DataArray -= Mu
            DataArray /= Sd
            (DStat,
             PVal)     = scipy.stats.kstest(DataArray,
                                            'norm')
            fw.write(str(DStat)+',')
            fw.write(str(PVal))
            fw.write('\n')
    fw.close()

    LogString = '  ** Performing normality tests on the ln-transformed concentrations ...'
    print LogString
    Log.Write(LogString+'\n')

    fw = open('NormalityCheckLnConcentrations.csv','w')
    fw.write('metabolite,KS test statistic, KS p-value\n')
    for Value in DataDict.itervalues():
        if(Value.GetMetaboliteName() and
           (not Value.GetMetaboliteName() in ExcludedMtbList)):
            fw.write(Value.GetMetaboliteConventionName()+',')

            DataArray = None
            if(Value.GetImputedDataArray()!=None):
                DataArray = scipy.real(scipy.array(Value.GetImputedDataArray()))
            else:
                DataArray  = scipy.real(scipy.array(Value.GetDataArray()))
            FilterArray = (DataArray!=0.0)
            DataArray   = scipy.compress(FilterArray,DataArray)
            DataArray   = scipy.log(DataArray)

            Mu         = DataArray.mean()
            Sd         = DataArray.std(ddof=1)
            DataArray -= Mu
            DataArray /= Sd
            (DStat,
             PVal)     = scipy.stats.kstest(DataArray,
                                            'norm')
            fw.write(str(DStat)+',')
            fw.write(str(PVal))
            fw.write('\n')
    fw.close()

    LogString = '  ** Performing normality tests on all ratios of ln-transformed concentrations ...'
    print LogString
    Log.Write(LogString+'\n')
    Keys = []
    for Key in DataDict.iterkeys():
        if(DataDict[Key].GetMetaboliteName() and
          (not DataDict[Key].GetMetaboliteName() in ExcludedMtbList)):
            Keys.append(Key)

    fw = open('NormalityCheckLnRatios.csv','w')
    fw.write('metaboliteratio,KS test statistic, KS p-value\n')
    for i in range(len(Keys)-1):
        KeyI = Keys[i]
        for j in range(i+1,len(Keys)):
            KeyJ = Keys[j]

            RatioName  = DataDict[KeyI].GetMetaboliteConventionName()
            RatioName += '_Over_'
            RatioName += DataDict[KeyJ].GetMetaboliteConventionName()
            fw.write(RatioName+',')

            boPass    = True
            DataArray = None
            if((DataDict[KeyI].GetImputedDataArray()!=None) and
               (DataDict[KeyJ].GetImputedDataArray()!=None)):
                ArrayI = scipy.real(scipy.array(DataDict[KeyI].GetImputedDataArray()))
                ArrayJ = scipy.real(scipy.array(DataDict[KeyJ].GetImputedDataArray()))

                FilterArray  = (ArrayI!=0.0)
                FilterArray *= (ArrayJ!=0.0)
                ArrayI       = scipy.compress(FilterArray,ArrayI)
                ArrayJ       = scipy.compress(FilterArray,ArrayJ)
                if(len(ArrayI)>0):
                    DataArray    = scipy.log(ArrayI)-scipy.log(ArrayJ)
                else:
                    boPass = False
            else:
                ArrayI = scipy.real(scipy.array(DataDict[KeyI].GetDataArray()))
                ArrayJ = scipy.real(scipy.array(DataDict[KeyJ].GetDataArray()))

                FilterArray  = (ArrayI!=0.0)
                FilterArray *= (ArrayJ!=0.0)
                ArrayI       = scipy.compress(FilterArray,ArrayI)
                ArrayJ       = scipy.compress(FilterArray,ArrayJ)
                if(len(ArrayI)>0):
                    DataArray    = scipy.log(ArrayI)-scipy.log(ArrayJ)
                else:
                    boPass = False

            if(boPass):
                Mu         = DataArray.mean()
                Sd         = DataArray.std(ddof=1)
                DataArray -= Mu
                DataArray /= Sd
                (DStat,
                 PVal)     = scipy.stats.kstest(DataArray,
                                                'norm')
                fw.write(str(DStat)+',')
                fw.write(str(PVal))
                fw.write('\n')
            else:
                fw.write('NA,')
                fw.write('NA')
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
    if(Arguments.ExcelReferenceDataFileName!=None):
        LogString = '**** Parsing Excel book \"'+Arguments.ExcelReferenceDataFileName+'\" containing the reference data ...'
        print LogString
        Log.Write(LogString+'\n')
        ReferenceDataDict = ExcelParser.ParseExcellBook(Arguments.ExcelReferenceDataFileName,
                                                        Arguments.ExcelReferenceSheetName,
                                                        Arguments.boRemoveCustomMetabolites,
                                                        Log)
    else:
        LogString = '**** No reference data file given: performing QC without reference data!'
        print LogString
        Log.Write(LogString+'\n')

    LogString = '**** Parsing analytical ranges from \"'+str(Arguments.AnalyticalRangesExcelFileName)+'\" ...'
    print LogString
    Log.Write(LogString+'\n')
    AnalRanges = BioCratesAnalyticalRanges.BioCratesAnalyticalRanges(Arguments.AnalyticalRangesExcelFileName)
    for Value in SampleDataDict.itervalues():
        Value.SetLODFromDocumentation(AnalRanges)
        Value.SetLLOQFromDocumentation(AnalRanges)
        Value.SetULOQFromDocumentation(AnalRanges)
    if(Arguments.ExcelReferenceDataFileName!=None):
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
    if((Arguments.boQCReference) and
       (Arguments.ExcelReferenceDataFileName!=None)):
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
    if(((Arguments.boQCReference) and
        (Arguments.ExcelReferenceDataFileName!=None)) or
       Arguments.boQCSample):
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
                                                     SampleDataDict,
                                                     Arguments.boImputeZeros)
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
