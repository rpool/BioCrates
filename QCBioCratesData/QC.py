import copy
import collections
import re
import scipy
import sys
import gzip
import argparse

import Logger
import MetabolomicsQCProtocol

class RMiceImpute:
    def __init__(self):
        self.RObjectFileName = None
        self.RReportFileName = None

    def SetRobjectFileName(self):
        self.RObjectFileName = 'MiceImputeRObject.robject'
    def GetRobjectFileName(self):
        if(self.RObjectFileName==None):
            self.SetRobjectFileName()
        return self.RObjectFileName
    def SetRReportFileName(self):
        self.RReportFileName = 'MiceImpute.report.gz'
    def GetRReportFileName(self):
        if(self.RReportFileName==None):
            self.SetRReportFileName()
        return self.RReportFileName
    def Impute(self,
               Arguments=argparse.Namespace,
               InputDataDict=dict,
               Log=Logger): # impute the missing values using the R-package mice
        import rpy2.robjects as robjects
        from rpy2.robjects.packages import importr
        import rpy2.rinterface as rinterface

        Mice    = None
        Imputed = None
        if(not Arguments.boReadRObject):
            LogString  = '** Writing rpy2 imputation output to \"'+Log.GetFileName()+'\" ...\n'
            LogString += '## START rpy2 ##'
            print LogString
            Log.Write(LogString+'\n')

            StdOutSav   = sys.stdout
            sys.stdout  = Log.GetFileHandle()
            DataFrameDict = {}
            for Key,Value in InputDataDict.iteritems():
                if(Value.GetMetaboliteName()!=None):
                    hlen = len(Value.GetDataArray())
                    DataFrameDict[str(Key)] = rinterface.baseenv['as.real'](rinterface.StrSexpVector(Value.GetDataArray()[0:hlen]))
                    RDataFrame  = robjects.DataFrame(DataFrameDict)
            Mice        = importr('mice')
            Imputed     = Mice.mice(RDataFrame)
    #===============================================================================
    # #        for testing purposes (takes a lot less time...)
    #        Imputed     = Mice.mice(Mice.nhanes)
    #===============================================================================
            sys.stdout  = StdOutSav
            LogString   = '## END rpy2 ##'
            print LogString
            Log.Write(LogString+'\n')

            LogString = '** Writing mice imputed R-object to \"'+self.GetRobjectFileName()+'\" ...'
            print LogString
            Log.Write(LogString+'\n')
            robjects.r.saveRDS(Imputed,file=self.GetRobjectFileName())
        else:
            LogString = '** Reading mice imputed R-object from \"'+self.GetRobjectFileName()+'\" ...'
            Mice      = importr('mice')
            Imputed   = robjects.r.readRDS(file=self.GetRobjectFileName())
            print LogString
            Log.Write(LogString+'\n')

        LogString = '** Writing mice imputed R-objects \'data\' and \'complete\' to \"'+self.GetRReportFileName()+'\" ...'
        print LogString
        Log.Write(LogString+'\n')
        fw = gzip.open(self.GetRReportFileName(),'wb')
        fw.write('## START TIMESTAMP\n')
        fw.write(str(Log.GetStartDate())+'\n')
        fw.write('## END TIMESTAMP\n')
        fw.write('## START DATA\n')
        DataNdArray = scipy.array(Imputed[Imputed.names.index('data')])
        i_imp       = Imputed.names.index('imp')
        imp         = Imputed[i_imp]
        DataDict    = {}
        for i in range(len(imp.names)):
            Name = imp.names[i]
            DataDict[Name] = DataNdArray[i]
        FormatString = '{0:>'+str(10)+'}'
        Line         = ''
        for i in range(len(imp.names)):
            Line += FormatString.format(imp.names[i])
        Line += '\n'
        fw.write(Line)
        for i in range(len(DataDict[imp.names[-1]])):
            Line = ''
            for j in range(len(imp.names)):
                Name  = imp.names[j]
                Entry = DataDict[Name][i]
                if(str(Entry)=='nan'):
                    Entry = 'NA'
                else:
                    Entry = str(round(Entry,3))
                Line += FormatString.format(Entry)
            Line += '\n'
            fw.write(Line)
#        fw.write(str(Imputed[Imputed.names.index('data')]))
        fw.write('## END DATA\n')
        fw.write('## START COMPLETE\n')
        Completed      = Mice.complete(Imputed,'repeated')
        CompletedNames = list(Completed.names)
        CompletedArray = scipy.array(Completed)
        i_imp          = Imputed.names.index('imp')
        imp            = Imputed[i_imp]
        DataDict       = {}
        for i in range(len(imp.names)):
            Name        = imp.names[i]
            AvArray     = None
            Counter     = 0
            for CName in CompletedNames:
                if(Name==CName.split('.')[0]):
                    Index = CompletedNames.index(CName)
                    if(AvArray==None):
                        AvArray = CompletedArray[Index]
                    else:
                        AvArray += CompletedArray[Index]
                    Counter += 1
            AvArray /= float(Counter)
            DataDict[Name] = AvArray
        FormatString = '{0:>'+str(10)+'}'
        Line         = ''
        for i in range(len(imp.names)):
            Line += FormatString.format(imp.names[i])
        Line += '\n'
        fw.write(Line)
        for i in range(len(DataDict[imp.names[-1]])):
            Line = ''
            for j in range(len(imp.names)):
                Name  = imp.names[j]
                Entry = DataDict[Name][i]
                if(str(Entry)=='nan'):
                    Entry = 'NA'
                else:
                    Entry = str(round(Entry,3))
                Line += FormatString.format(Entry)
            Line += '\n'
            fw.write(Line)
        fw.write('## END COMPLETE\n')
        fw.close()

        return

class QualityControl:
    def __init__(self):
        self.MetaboliteExclusionList     = None
        self.SampleExclusionList         = None
        self.MapMetaboliteListToDataDict = None
        self.MetabolitesWitMissings      = None
        return
    def SetMetaboliteExclusionList(self,
                                   MetaboliteContainers=list):
        self.MetaboliteExclusionList = []
        for i in range(len(MetaboliteContainers)):
            if(MetaboliteContainers[i].GetExcludeMetabolite()):
                self.MetaboliteExclusionList.append(i)
        self.MetaboliteExclusionList.sort()
        self.MetaboliteExclusionList.reverse()
        return
    def GetMetaboliteExclusionList(self,
                                   MetaboliteContainers=list):
        self.SetMetaboliteExclusionList(MetaboliteContainers)
        return self.MetaboliteExclusionList
    def SetSampleExclusionList(self,
                               MetaboliteContainers=list):
        self.SampleExclusionList = []
        for i in range(len(MetaboliteContainers)):
            self.SampleExclusionList.append(MetaboliteContainers[i].GetSampleExclusionList())
        if(len(MetaboliteContainers)==len(self.SampleExclusionList)):
            TmpSampleExclusionList = []
            for Entry in self.SampleExclusionList:
                if(not (Entry in TmpSampleExclusionList)):
                    TmpSampleExclusionList.append(Entry)
            if(len(TmpSampleExclusionList)==1):
                self.SampleExclusionList = TmpSampleExclusionList[0]
            if(self.SampleExclusionList==None):
                self.SampleExclusionList = []
        return
    def GetSampleExclusionList(self,
                               MetaboliteContainers=list):
        self.SetSampleExclusionList(MetaboliteContainers)
        return self.SampleExclusionList
    def DetectMetabolitesWitMissings(self,
                                     MetaboliteContainers=list):
        self.MetabolitesWitMissings = []
        for i in range(len(MetaboliteContainers)):
            if(MetaboliteContainers[i].HasMissingValues()):
                self.MetabolitesWitMissings.append(i)
        return self.MetabolitesWitMissings


def QCReference(DataDict={},                        # should contain keys of type int and values of type DataContainer
                Log=Logger):
    QCProtocol = MetabolomicsQCProtocol.MetabolomicsQCProtocol()

    MetaboliteContainers = []
    for Key,Value in DataDict.iteritems():
        if(Value.GetDataName()=='Plate bar code'):
            PlateBarCodeContainer = copy.deepcopy(Value)
        elif(Value.GetMetaboliteName()!=None):
            MetaboliteContainers.append(copy.deepcopy(Value))

    for i in range(len(MetaboliteContainers)):
        MetaboliteContainers[i].SetExcludeMetabolite(False)

#   ENGAGE QC V2.0 1a)
    for i in range(len(MetaboliteContainers)):
        PlateBarCodeDataArray        = PlateBarCodeContainer.GetDataArray()
        MetaboliteDataArray          = MetaboliteContainers[i].GetDataArray()
        MetaboliteName               = MetaboliteContainers[i].GetMetaboliteName()
        ConcentrationPerPlateBarCode = {}
        for j in range(len(PlateBarCodeDataArray)):
            PlateBarCode = PlateBarCodeDataArray[j]
            if(not ConcentrationPerPlateBarCode.has_key(PlateBarCode)):
                ConcentrationPerPlateBarCode[PlateBarCode] = []
            ConcentrationPerPlateBarCode[PlateBarCode].append(MetaboliteDataArray[j])
        TotalCV = 0.0
        NTot    = 0
        for Key,Value in ConcentrationPerPlateBarCode.iteritems():
            Std      = scipy.array(Value).std(ddof=1)
            Mean     = scipy.array(Value).mean()
            CV       = Std/Mean
            print '****** ',MetaboliteName,len(Value)
            if(str(CV)==str(scipy.nan) or
               str(CV)==str(scipy.inf)): # 0/0 or x/0 (x>0) division
                LogString  = '** Zero division for CV of \"'+MetaboliteName+'\"/'+str(Key)+' Metabolite/PlateBarCode combination ...'
                print LogString
                Log.Write(LogString+'\n')
                if(not MetaboliteContainers[i].GetExcludeMetabolite()):
                    LogString = '** Setting exclusion status of metabolite \"'+MetaboliteName+'\" to True in ReferenceDataDict ...'
                    print LogString
                    Log.Write(LogString+'\n')
                    MetaboliteContainers[i].SetExcludeMetabolite(True)
            else:
                n        = len(Value)
                TotalCV += CV
                NTot    += 1
        if(not MetaboliteContainers[i].GetExcludeMetabolite()):
            MeanTotalCV = 100.0*TotalCV/float(NTot)
            if(MeanTotalCV>QCProtocol.GetRefererenceMeanCVOfAllPlatesThreshold()):
                LogString  = '** Mean(CV of all plates) > '+str(QCProtocol.GetRefererenceMeanCVOfAllPlatesThreshold())+'% for Metabolite '+MetaboliteName+' ...\n'
                LogString += '** Mean(CV of all plates) = '+str(round(MeanTotalCV,3))+'% ...\n'
                LogString += '** Setting exclusion status of metabolite \"'+MetaboliteName+'\" to True in ReferenceDataDict ...'
                print LogString
                Log.Write(LogString+'\n')
                MetaboliteContainers[i].SetExcludeMetabolite(True)

    return MetaboliteContainers # remember that these datacontainers are copies,
                                # they do not affect the input object objects

def QCSample(DataDict={},                        # should contain keys of type int and values of type DataContainer
             ReferenceQCMetaboliteContainers=[], # should contain entries of type DataContainer
             DLHandling=str,                     # this determines how the detection limit properties of the metabolilte are handled
             Log=Logger):
    QCProtocol = MetabolomicsQCProtocol.MetabolomicsQCProtocol()

    MetaboliteContainers                  = []
    MapMetaboliteContainersToDataDictKeys = []
    for Key,Value in DataDict.iteritems():
        if(Value.GetMetaboliteName()!=None):
            MetaboliteContainers.append(copy.deepcopy(Value))
            MapMetaboliteContainersToDataDictKeys.append(Key)

    for i in range(len(MetaboliteContainers)):
        MetaboliteContainers[i].SetExcludeMetabolite(False)
        MetaboliteContainers[i].SetMeanMetaboliteConcentrationAndStd()

#   ENGAGE QC V2.0 1a) Check for concentrations below LOD
    if(re.search('LOD',DLHandling)):
        for i in range(len(MetaboliteContainers)):
            MetaboliteDataArray = MetaboliteContainers[i].GetDataArray()
            MetaboliteName      = MetaboliteContainers[i].GetMetaboliteName()
            MetaboliteLODValue  = MetaboliteContainers[i].GetLODFromDocumentation()
            MetaboliteLLOQValue = MetaboliteContainers[i].GetLLOQFromDocumentation()
            boValuesBelowLOD    = False
            MissingValueRate    = 0
            for j in range(len(MetaboliteDataArray)):
                if(MetaboliteDataArray[j]<MetaboliteLODValue):
                    boValuesBelowLOD       = True
                    MetaboliteDataArray[j] = MetaboliteContainers[i].GetMissingDataIdentifier()
                    MissingValueRate      += 1
            if(boValuesBelowLOD):
                LogString  = '** The metabolite data of '+MetaboliteName+' contain values below the LOD ...\n'
                LogString += '** These values were set to \"'+MetaboliteContainers[i].GetMissingDataIdentifier()+'\" in the SampleDataDict ...\n'
                LogString += '** Missing value rate = '+str(round(100.0*float(MissingValueRate)/float(len(MetaboliteDataArray)),3))+'% ...'
                if(100.0*float(MissingValueRate)/float(len(MetaboliteDataArray))>QCProtocol.GetMetaboliteMissingValueRateThreshold()):
                    LogString += '\n** Missing value rate > '+str(QCProtocol.GetMetaboliteMissingValueRateThreshold())+'% => setting exclusion status of metabolite to True in SampleDataDict ...'
                    MetaboliteContainers[i].SetExcludeMetabolite(True)
                print LogString
                Log.Write(LogString+'\n')
            if((DLHandling=='LOD-LLOQ') and
               (MetaboliteLLOQValue!=None) and
               (not MetaboliteContainers[i].GetExcludeMetabolite())): # Exclude metabolite if mean(concentration) < LLOQ
                MeanConc,\
                MeanStd    = MetaboliteContainers[i].GetMeanMetaboliteConcentrationAndStdExcludingMissingValues()
                if(MeanConc<MetaboliteLLOQValue):
                    LogString  = '** The mean metabolite concentration of '+MetaboliteName+' is below the LLOQ ...\n'
                    LogString += '** These values are set to \"'+MetaboliteContainers[i].GetMissingDataIdentifier()+'\" in the SampleDataDict ...\n'
                    LogString += '** Setting exclusion status of metabolite to True in SampleDataDict ...'
                    print LogString
                    Log.Write(LogString+'\n')
                    MetaboliteContainers[i].SetExcludeMetabolite(True)
#   ENGAGE QC V2.0 1b) Check for outliers and impute missing data
    for i in range(len(MetaboliteContainers)):
        MetaboliteDataArray = MetaboliteContainers[i].GetDataArray()
        MetaboliteName      = MetaboliteContainers[i].GetMetaboliteName()
        MeanConcentration,\
        StdConcentration    = MetaboliteContainers[i].GetMeanMetaboliteConcentrationAndStdExcludingMissingValues()
        OutlierIndexList    = []
        for j in range(len(MetaboliteDataArray)):
            Concentration = MetaboliteDataArray[j]
            if(Concentration!=MetaboliteContainers[i].GetMissingDataIdentifier()):
                DeviationInSigma = (Concentration - MeanConcentration)/StdConcentration
                if(abs(DeviationInSigma)>QCProtocol.GetMetaboliteOutlyingSampleThreshold()):
                    LogString  = '** Detected a metabolite concentration outlier for \"'+MetaboliteName+'\", (index '+str(j)+')...\n'
                    LogString += '** Setting this entry to \"'+MetaboliteContainers[i].GetMissingDataIdentifier()+'\" ...\n'
                    LogString += '** Appending this entry to the OutlierIndexList ...'
                    print LogString
                    Log.Write(LogString+'\n')
                    MetaboliteContainers[i].AlterDataArrayEntry(j,
                                                                MetaboliteContainers[i].GetMissingDataIdentifier())
                    OutlierIndexList.append(j)
        MetaboliteContainers[i].SetOutlierIndexList(OutlierIndexList)
#   Detect outliers per sample
    DuplicateList = []
    for i in range(len(MetaboliteContainers)):
        ReferenceboExcludeMetabolite = None
        if(len(ReferenceQCMetaboliteContainers)==0):
            ReferenceboExcludeMetabolite = False
        else:
            ReferenceboExcludeMetabolite = ReferenceQCMetaboliteContainers[i].GetExcludeMetabolite()
        if(not (MetaboliteContainers[i].GetExcludeMetabolite() or
                ReferenceboExcludeMetabolite)):
            OutLierList = MetaboliteContainers[i].GetOutlierIndexList()
            DuplicateList.extend(OutLierList)
    DuplicateDict    = collections.Counter(DuplicateList)
    SetToMissingList = []
    for Key, Value in DuplicateDict.iteritems(): # Key=DataArrayIndex;Val=ListOfMetaboliteContainers
        if(Value<=3):
            SetToMissingList.append(Key)
    for Index in SetToMissingList:
        del DuplicateDict[Index]

    ExcludeSampleList = []
    for Key,Value in DuplicateDict.iteritems(): # Key=DataArrayIndex;Val=ListOfMetaboliteContainers
        DataArrays = []
        for i in range(len(MetaboliteContainers)): # get DataArrays containing the FULL raw data
            if(Key in MetaboliteContainers[i].GetOutlierIndexList()):
                DataKey   = MapMetaboliteContainersToDataDictKeys[i]
                DataArrays.append(scipy.array(DataDict[DataKey].GetDataArray()))
        DataArrays     = scipy.array(DataArrays)
        R2NdArray      = scipy.power(scipy.corrcoef(DataArrays),2.0) # R^2 is is the proportion of variance in Y that can be accounted for by knowing X, and vice-versa
        OutlierCounter = 0
        for i in range(len(R2NdArray)):
            boOutlier = True
            for j in range(len(R2NdArray[i])):
                if((i!=j) and
                   (100.0*R2NdArray[i,j]>=QCProtocol.GetOutlierCorrelationPercentage())):
                    boOutlier = False
            if(boOutlier):
                OutlierCounter +=1
        if(OutlierCounter>3):
            LogString  = '** The number of independent outliers for sample '+str(Key)+' exceeds 3 ...\n'
            LogString += '** Such an outlier is defined by an \'outlying metabolite observation that has a correlation of less than\n'
            LogString += '** '+str(QCProtocol.GetOutlierCorrelationPercentage())+'% with all other outlying metabolite concentrations of that sample ...\n'
            LogString += '** Sample '+str(Key)+' will be appended to the sample exclusion list ...'
            print LogString
            Log.Write(LogString+'\n')
            ExcludeSampleList.append(Key)
    for i in range(len(MetaboliteContainers)):
        MetaboliteContainers[i].ExtendToSampleExclusionList(ExcludeSampleList)
    for Index in SetToMissingList:
        for i in range(len(MetaboliteContainers)):
            if(Index in MetaboliteContainers[i].GetOutlierIndexList()):
                MetaboliteName = MetaboliteContainers[i].GetMetaboliteName()
                LogString  = '** Setting the DataValue of '+MetaboliteName+ '(index: '+str(Index)+') to missing ...'
                print LogString
                Log.Write(LogString+'\n')
                MetaboliteContainers[i].AlterDataArrayEntry(Index,
                                                            MetaboliteContainers[i].GetMissingDataIdentifier())

    return MetaboliteContainers # remember that these datacontainers are copies,
                                # they do not affect the input object objects
def ProcessDataDict(ReferenceQCMetaboliteContainers=list,
                    SampleQCMetaboliteContainers=list,
                    SampleDataDict=dict,
                    boImputeZeros=bool):
    QCntrl                                = QualityControl()
    MetaboliteExclusionList               = []
    SampleExclusionList                   = []
    MapMetaboliteContainersToDataDictKeys = []
    MetaboliteNameExclusionList           = []
    SampleIDExclusionList                 = []

    for Key,Value in SampleDataDict.iteritems():
            if(Value.GetMetaboliteName()!=None):
                MapMetaboliteContainersToDataDictKeys.append(Key)
    if(len(ReferenceQCMetaboliteContainers)>0):
        MetaboliteExclusionList.extend(QCntrl.GetMetaboliteExclusionList(ReferenceQCMetaboliteContainers))
    MetaboliteExclusionList.extend(QCntrl.GetMetaboliteExclusionList(SampleQCMetaboliteContainers))
    MetaboliteExclusionList = list(set(MetaboliteExclusionList))
    MetaboliteExclusionList.sort()
    for i in MetaboliteExclusionList:
        MetaboliteNameExclusionList.append(SampleDataDict[MapMetaboliteContainersToDataDictKeys[i]].GetMetaboliteName())
    MetaboliteExclusionList.reverse()

    if(len(ReferenceQCMetaboliteContainers)>0):
        SampleExclusionList.extend(QCntrl.GetSampleExclusionList(ReferenceQCMetaboliteContainers))
    SampleExclusionList.extend(QCntrl.GetSampleExclusionList(SampleQCMetaboliteContainers))
    SampleExclusionList = list(set(SampleExclusionList))
    SampleExclusionList.sort()
    SampleIDKeyInDataDict = None
    for Key,Value in SampleDataDict.iteritems():
        if(Value.GetDataName()=='Sample Identification'):
            SampleIDKeyInDataDict = Key
            break
    for i in SampleExclusionList:
        SampleIDExclusionList.append(SampleDataDict[SampleIDKeyInDataDict].GetDataArray()[i])
    SampleExclusionList.reverse()

    ProcessedDataDict = copy.deepcopy(SampleDataDict)
    for i in range(len(SampleQCMetaboliteContainers)):
        for j in SampleQCMetaboliteContainers[i].GetOutlierIndexList():
            ProcessedDataDict[MapMetaboliteContainersToDataDictKeys[i]].\
            AlterDataArrayEntry(j,
                                ProcessedDataDict[MapMetaboliteContainersToDataDictKeys[i]].GetMissingDataIdentifier())
        if(boImputeZeros):
            for j in range(len(SampleQCMetaboliteContainers[i].GetDataArray())):
                if((type(SampleQCMetaboliteContainers[i].GetDataArray()[j])==float) and
                   (SampleQCMetaboliteContainers[i].GetDataArray()[j]<=0.0)):
                    ProcessedDataDict[MapMetaboliteContainersToDataDictKeys[i]].\
                    AlterDataArrayEntry(j,
                                        ProcessedDataDict[MapMetaboliteContainersToDataDictKeys[i]].GetMissingDataIdentifier())
    for i in MetaboliteExclusionList:
        del ProcessedDataDict[MapMetaboliteContainersToDataDictKeys[i]]
    for Key,Value in ProcessedDataDict.iteritems():
        for i in SampleExclusionList:
            Value.DeleteIndexFromDataArray(i)

    return ProcessedDataDict,\
           MetaboliteNameExclusionList,\
           SampleIDExclusionList

def QCFinal(InputDataDict=dict,
            Log=Logger):
    QCProtocol                  = MetabolomicsQCProtocol.MetabolomicsQCProtocol()
    OutputDataDict              = copy.deepcopy(InputDataDict)
    ExcludeMetaboliteList       = []
    MetaboliteNameExclusionList = []
    for Key,Value in OutputDataDict.iteritems():
        if(Value.GetMetaboliteName()!=None):
            Counter = 0
            for Entry in Value.GetDataArray():
                if(Entry==Value.GetMissingDataIdentifier()):
                    Counter+=1
            MissingRate = float(Counter/float(len(Value.GetDataArray())))
            if(MissingRate*100.0>QCProtocol.GetMetaboliteMissingValueRateThreshold()):
                LogString  = '** Missing value rate for metabolite \"'+Value.GetMetaboliteName()+'\" is too high ...\n'
                LogString += '** ('+str(round(MissingRate*100.0,3))+'% > '+str(QCProtocol.GetMetaboliteMissingValueRateThreshold())+\
                            '%) ...\n'
                LogString += '** Excluding this metabolite from analysis ...'
                print LogString
                Log.Write(LogString+'\n')
                ExcludeMetaboliteList.append(Key)
    ExcludeMetaboliteList.sort()
    for i in ExcludeMetaboliteList:
        MetaboliteNameExclusionList.append(InputDataDict[i].GetMetaboliteName())
    ExcludeMetaboliteList.reverse()
    for i in ExcludeMetaboliteList:
        del OutputDataDict[i]

    return OutputDataDict,\
           MetaboliteNameExclusionList