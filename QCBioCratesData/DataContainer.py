import sys
import scipy
import copy
import re

import BioCratesAnalyticalRanges

class DataContainer:
    """
    Class containing the actual data.
    """
    def __init__(self):
        self.DataName                                          = None
        self.MetaboliteClass                                   = None
        self.MetaboliteName                                    = None
        self.MetaboliteConventionName                          = None
        self.MetaboliteIndex                                   = None
        self.LODInfoArray                                      = None
        self.LODValueArray                                     = None
        self.DataArray                                         = None
        self.QCedDataArray                                     = None
        self.ImputedDataArray                                  = None
        self.boExcludeMetabolite                               = None # only if mean(CV of all plates) > 25%;
                                                                      # if type==None: this is not a metabolite
        self.LODFromDocumentation                              = None
        self.LLOQFromDocumentation                             = None
        self.ULOQFromDocumentation                             = None
        self.MeanMetaboliteConcentration                       = None
        self.StdMetaboliteConcentraion                         = None
        self.MeanMetaboliteConcentrationExcludingMissingValues = None
        self.StdMetaboliteConcentraionExlcudingMissingValues   = None
        self.MissingDataIdentifier                             = None
        self.SetMissingDataIdentifier()
        self.OutlierIndexList                                  = None # A list of outliers in reverse sorted order
        self.SampleExclusionList                               = None # A list of samples (persons that should be excluded)

    def SetDataName(self,
                    String=str):
        self.DataName = String
        return
    def GetDataName(self):
        return self.DataName
    def SetMetaboliteClass(self,
                           String=str):
        self.MetaboliteClass = String
        return
    def GetMetaboliteClass(self):
        return self.MetaboliteClass
    def SetMetaboliteName(self,
                           String=str):
        self.MetaboliteName = String
        return
    def GetMetaboliteName(self):
        return self.MetaboliteName
    def SetMetaboliteConventionName(self):
        self.MetaboliteConventionName = re.sub('[(,),\-,\ ,:]','.',self.GetMetaboliteName())
        self.MetaboliteConventionName = re.sub('\.\.','.',self.GetMetaboliteConventionName())
        if(self.GetMetaboliteConventionName()[-1]=='.'):
            self.MetaboliteConventionName = self.GetMetaboliteConventionName()[:-1]
        return
    def GetMetaboliteConventionName(self):
        if(self.MetaboliteConventionName==None):
            self.SetMetaboliteConventionName()
        return self.MetaboliteConventionName
    def InitLODInfoArray(self):
        if(self.LODInfoArray==None):
            self.LODInfoArray = []
        return
    def AppendToLODInfoArray(self,
                             String=str):
        self.InitLODInfoArray()
        self.LODInfoArray.append(String)
        return
    def GetLODInfoArray(self):
        return self.LODInfoArray
    def InitLODValueArray(self):
        if(self.LODValueArray==None):
            self.LODValueArray = []
        return
    def AppendToLODValueArray(self,
                              Value=float):
        self.InitLODValueArray()
        self.LODValueArray.append(Value)
        return
    def GetLODValueArray(self):
        return self.LODValueArray
    def InitDataArray(self):
        if(self.DataArray==None):
            self.DataArray = []
        return
    def AppendToDataArray(self,
                          Value=str):
        self.InitDataArray()
        self.DataArray.append(Value)
        return
    def GetDataArray(self):
        if(self.DataArray==None):
            print '!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            print '! ERROR DataArray not set !'
            print '!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            print 'EXITING ...'
            sys.exit(1)
        return self.DataArray
    def SetQCedDataArray(self,
                         Array=list):
        self.QCedDataArray = copy.deepcopy(Array)
        return
    def GetQCedDataArray(self):
        return self.QCedDataArray
    def SetImputedDataArray(self,
                            Array=list):
        self.ImputedDataArray = copy.deepcopy(Array)
        return
    def GetImputedDataArray(self):
        return self.ImputedDataArray
    def DeleteIndexFromDataArray(self,
                                 index=int):
        if(self.DataArray!=None):
            del self.DataArray[index]
        return
    def SetExcludeMetabolite(self,
                             Value=bool):
        self.boExcludeMetabolite = Value
        return
    def GetExcludeMetabolite(self):
        return self.boExcludeMetabolite
    def SetLODFromDocumentation(self,
                                AnalytRanges=BioCratesAnalyticalRanges):
        if(self.GetMetaboliteName()!=None):
            self.LODFromDocumentation = AnalytRanges.GetLODValueDict()[self.GetMetaboliteName()]
        return
    def GetLODFromDocumentation(self):
        return self.LODFromDocumentation
    def SetLLOQFromDocumentation(self,
                                 AnalytRanges=BioCratesAnalyticalRanges):
        if(self.GetMetaboliteName()!=None):
            self.LLOQFromDocumentation = AnalytRanges.GetLLOQValueDict()[self.GetMetaboliteName()]
        return
    def GetLLOQFromDocumentation(self):
        return self.LLOQFromDocumentation
    def SetULOQFromDocumentation(self,
                                 AnalytRanges=BioCratesAnalyticalRanges):
        if(self.GetMetaboliteName()!=None):
            self.ULOQFromDocumentation = AnalytRanges.GetULOQValueDict()[self.GetMetaboliteName()]
        return
    def GetULOQFromDocumentation(self):
        return self.ULOQFromDocumentation
    def SetMeanMetaboliteConcentrationAndStd(self):
        if(self.MetaboliteName!=None):
            DataArray = scipy.array(self.DataArray)
            self.MeanMetaboliteConcentration = DataArray.mean()
            self.StdMetaboliteConcentraion   = DataArray.std()
        return
    def GetMeanMetaboliteConcentrationAndStd(self):
        if(self.MeanMetaboliteConcentration==None):
            self.SetMeanMetaboliteConcentrationAndStd()
        return self.MeanMetaboliteConcentration,\
               self.StdMetaboliteConcentraion
    def SetMeanMetaboliteConcentrationAndStdExcludingMissingValues(self):
        boFloats  = scipy.array(list(type(Entry)==float for Entry in self.DataArray))
        TmpDataArray = scipy.compress(boFloats,scipy.array(self.DataArray))
        TmpDataArray = scipy.array(TmpDataArray,dtype=float)
        self.MeanMetaboliteConcentrationExcludingMissingValues = TmpDataArray.mean()
        self.StdMetaboliteConcentrationExlcudingMissingValues  = TmpDataArray.std()
        return
    def GetMeanMetaboliteConcentrationAndStdExcludingMissingValues(self):
        if(self.MetaboliteName!=None):
            self.SetMeanMetaboliteConcentrationAndStdExcludingMissingValues()
        return self.MeanMetaboliteConcentrationExcludingMissingValues,\
               self.StdMetaboliteConcentrationExlcudingMissingValues
    def SetMissingDataIdentifier(self):
        self.MissingDataIdentifier = 'NA'
        return
    def GetMissingDataIdentifier(self):
        if(self.MissingDataIdentifier==None):
            self.SetMissingDataIdentifier()
        return self.MissingDataIdentifier
    def SetOutlierIndexList(self,
                            OutlierIndexList=list):
        self.OutlierIndexList = []
        self.OutlierIndexList.extend(OutlierIndexList)
        self.OutlierIndexList.sort()
        self.OutlierIndexList.reverse()
        return
    def GetOutlierIndexList(self):
        if(self.OutlierIndexList==None):
            print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            print '! ERROR OutlierIndexList not set !'
            print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            print 'EXITING ...'
            sys.exit(1)
        return self.OutlierIndexList
    def AlterDataArrayEntry(self,
                            Index,
                            Value):
        if(self.DataArray==None):
            print '!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            print '! ERROR DataArray not set !'
            print '!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            print 'EXITING ...'
            sys.exit(1)
        self.DataArray[Index] = Value
        return
    def InitSampleExclusionList(self):
        self.SampleExclusionList = []
        return
    def AppendToSampleExclusionList(self,
                                    Index):
        try:
            (type(Index)==int)
        except:
            raise TypeError
        if(self.SampleExclusionList==None):
            self.InitSampleExclusionList()
        self.SampleExclusionList.append(Index)
        self.SampleExclusionList.sort()
        self.SampleExclusionList.reverse() # should be sorted and reversed
        return
    def ExtendToSampleExclusionList(self,
                                    List):
        try:
            (type(List)==list)
        except:
            raise TypeError
        if(self.SampleExclusionList==None):
            self.InitSampleExclusionList()
        self.SampleExclusionList.extend(List)
        self.SampleExclusionList.sort()
        self.SampleExclusionList.reverse() # should be sorted and reversed
        return
    def GetSampleExclusionList(self):
        return self.SampleExclusionList
    def HasMissingValues(self):
        return (self.GetMissingDataIdentifier() in self.GetDataArray())
    def SetMetaboliteIndex(self,
                           Index=int):
        self.MetaboliteIndex = Index
        return
    def GetMetaboliteIndex(self):
        return self.MetaboliteIndex

