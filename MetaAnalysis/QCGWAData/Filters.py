import lxml.etree
import scipy
import os
import re
import collections

import FilterFunction
import DataContainer

class Filters:
    def __init__(self):
        self.InitialArrayLength    = 0
        self.NDeletedDuplicateSNPs = 0
        self.MaxNDuplicateSNPs     = 0
        self.boDuplicateSNPWarning = False
        self.FilterFunctionDict    = None
        self.FilterReportDict      = None
        self.FilterReportPath      = None

        self.InitFilterReportDict()
        return

    def InitFilterReportDict(self):
        self.FilterReportDict = {}
        return

    def GetFilterReportDictDict(self):
        return self.FilterReportDict

    def SetFilterReportPath(self,
                            Value=str):
        Cwd                   = os.getcwd()
        self.FilterReportPath = os.path.join(Cwd,Value)
        if(not os.path.isdir(self.FilterReportPath)):
            os.mkdir(self.FilterReportPath)
        return

    def GetFilterReportPath(self):
        if(self.FilterReportPath==None):
            self.SetFilterReportPath('Filters')
        return self.FilterReportPath

    def WriteFilterReport(self,
                          FileName=str,
                          Tag=str):
        FilePath = os.path.join(self.GetFilterReportPath(),FileName)

        FH = open(FilePath,'w')
        for Key in self.GetFilterReportDictDict()[Tag].iterkeys():
            String = '**'+self.GetFilterReportDictDict()[Tag][Key]
            String = re.sub('\n','\n**',String)
            FH.write(String+'\n')
        FH.close()
        return

    def WriteCustomFilterReport(self,
                                FileName=str,
                                String=str):
        FilePath = os.path.join(self.GetFilterReportPath(),FileName)

        FH = open(FilePath,'w')
        FH.write(String+'\n')
        FH.close()
        return

    def InitFilterReportDictDict(self,
                                 Tag=str):
        self.FilterReportDict[Tag] = {}
        return

    def SetFilterReportDictDict(self,
                                ParentTag=str,
                                ChildTag=str,
                                Value=str):
        self.FilterReportDict[ParentTag][ChildTag] = Value
        return

    def FilterSEs(self,
                  XmlObj=lxml.etree._ElementTree,
                  DCs=DataContainer.DataContainers,
                  ColumnTag=str,
                  boDryRun=False):

        FiltersTag = None
        if(boDryRun):
            FiltersTag = 'DryRunFilters'
        else:
            FiltersTag = 'Filters'

        FilterTags  = XmlObj.getroot().find('MtbGWAColumns').find(ColumnTag).find(FiltersTag).text.split(',')

        self.InitFilterReportDictDict(ColumnTag)

        FilterArray = None
        for Tag in FilterTags:
            Operator     = XmlObj.getroot().find('QCFilters').find(Tag).find('Operator').text
            CompareValue = XmlObj.getroot().find('QCFilters').find(Tag).find('Compare').text
            ValueType    = XmlObj.getroot().find('QCFilters').find(Tag).find('CompareType').text

            FFunction = FilterFunction.FilterFunction(OperatorString=Operator,
                                                      CompareString=CompareValue,
                                                      CompareType=ValueType)

            InitLength  = len(DCs.DataContainers[ColumnTag].GetDataArray())
            FinalLength = None
            FilterArray = FFunction.Run(DataArray=DCs.DataContainers[ColumnTag].GetDataArray())

            if(not boDryRun):
                for Key in DCs.DataContainers.iterkeys():
                    DataArray = scipy.compress(FilterArray,
                                               DCs.DataContainers[Key].GetDataArray())
                    DCs.DataContainers[Key].ReplaceDataArray(DataArray)
                FinalLength = len(DCs.DataContainers[ColumnTag].GetDataArray())
            else:
                CounterDict = collections.defaultdict(int)
                for Entry in FilterArray:
                    CounterDict[Entry] += 1
                Difference  = CounterDict[False]
                FinalLength = InitLength-Difference


            self.SetFilterReportDictDict(ParentTag=ColumnTag,
                                         ChildTag=Tag,
                                         Value='Column Tag   = '+ColumnTag+'\n'+\
                                               'Filter Tag   = '+Tag+'\n'+\
                                               'Start Length = '+str(InitLength)+'\n'+\
                                               'Final Length = '+str(FinalLength)+'\n'+\
                                               'Difference   = '+str(InitLength-FinalLength))

        return DCs,\
               FilterTags

    def FilterBetas(self,
                    XmlObj=lxml.etree._ElementTree,
                    DCs=DataContainer.DataContainers,
                    ColumnTag=str,
                    boDryRun=False):

        FiltersTag = None
        if(boDryRun):
            FiltersTag = 'DryRunFilters'
        else:
            FiltersTag = 'Filters'

        FilterTags  = XmlObj.getroot().find('MtbGWAColumns').find(ColumnTag).find(FiltersTag).text.split(',')

        self.InitFilterReportDictDict(ColumnTag)

        FilterArray = None
        for Tag in FilterTags:
            Operator     = XmlObj.getroot().find('QCFilters').find(Tag).find('Operator').text
            CompareValue = XmlObj.getroot().find('QCFilters').find(Tag).find('Compare').text
            ValueType    = XmlObj.getroot().find('QCFilters').find(Tag).find('CompareType').text

            FFunction = FilterFunction.FilterFunction(OperatorString=Operator,
                                                      CompareString=CompareValue,
                                                      CompareType=ValueType)

            InitLength  = len(DCs.DataContainers[ColumnTag].GetDataArray())
            FinalLength = None
            FilterArray = FFunction.Run(DataArray=DCs.DataContainers[ColumnTag].GetDataArray())

            if(not boDryRun):
                for Key in DCs.DataContainers.iterkeys():
                    DataArray = scipy.compress(FilterArray,
                                               DCs.DataContainers[Key].GetDataArray())
                    DCs.DataContainers[Key].ReplaceDataArray(DataArray)
                FinalLength = len(DCs.DataContainers[ColumnTag].GetDataArray())
            else:
                CounterDict = collections.defaultdict(int)
                for Entry in FilterArray:
                    CounterDict[Entry] += 1
                Difference  = CounterDict[False]
                FinalLength = InitLength-Difference

            self.SetFilterReportDictDict(ParentTag=ColumnTag,
                                         ChildTag=Tag,
                                         Value='Column Tag   = '+ColumnTag+'\n'+\
                                               'Filter Tag   = '+Tag+'\n'+\
                                               'Start Length = '+str(InitLength)+'\n'+\
                                               'Final Length = '+str(FinalLength)+'\n'+\
                                               'Difference   = '+str(InitLength-FinalLength))

        return DCs,\
               FilterTags

    def FilterChrs(self,
                   XmlObj=lxml.etree._ElementTree,
                   DCs=DataContainer.DataContainers,
                   ColumnTag=str):

        FilterTags  = XmlObj.getroot().find('MtbGWAColumns').find(ColumnTag).find('Filters').text.split(',')

        self.InitFilterReportDictDict(ColumnTag)

        FilterArray = None
        for Tag in FilterTags:
            Operator      = XmlObj.getroot().find('QCFilters').find(Tag).find('Operator').text
            CompareValues = XmlObj.getroot().find('QCFilters').find(Tag).find('Compare').text.split(',')
            ValueType     = XmlObj.getroot().find('QCFilters').find(Tag).find('CompareType').text

            FilterArray = None
            for i in range(len(CompareValues)):
                CValue    = CompareValues[i]
                FFunction = FilterFunction.FilterFunction(OperatorString=Operator,
                                                          CompareString=CValue,
                                                          CompareType=ValueType)
                if(i==0):
                    FilterArray = FFunction.Run(DataArray=DCs.DataContainers[ColumnTag].GetDataArray())
                else:
                    FilterArray = (FilterArray | FFunction.Run(DataArray=DCs.DataContainers[ColumnTag].GetDataArray()))

            InitLength  = len(DCs.DataContainers[ColumnTag].GetDataArray())
            for Key in DCs.DataContainers.iterkeys():
                DataArray = scipy.compress(FilterArray,
                                           DCs.DataContainers[Key].GetDataArray())
                DCs.DataContainers[Key].ReplaceDataArray(DataArray)
            FinalLength = len(DCs.DataContainers[ColumnTag].GetDataArray())

            self.SetFilterReportDictDict(ParentTag=ColumnTag,
                                         ChildTag=Tag,
                                         Value='Column Tag   = '+ColumnTag+'\n'+\
                                               'Filter Tag   = '+Tag+'\n'+\
                                               'Start Length = '+str(InitLength)+'\n'+\
                                               'Final Length = '+str(FinalLength)+'\n'+\
                                               'Difference   = '+str(InitLength-FinalLength))

        return DCs,\
               FilterTags

    def FilterNTotals(self,
                      XmlObj=lxml.etree._ElementTree,
                      DCs=DataContainer.DataContainers,
                      ColumnTag=str,
                      boDryRun=False):

        FiltersTag = None
        if(boDryRun):
            FiltersTag = 'DryRunFilters'
        else:
            FiltersTag = 'Filters'

        FilterTags  = XmlObj.getroot().find('MtbGWAColumns').find(ColumnTag).find(FiltersTag).text.split(',')

        self.InitFilterReportDictDict(ColumnTag)

        FilterArray = None
        for Tag in FilterTags:
            Operator     = XmlObj.getroot().find('QCFilters').find(Tag).find('Operator').text
            CompareValue = XmlObj.getroot().find('QCFilters').find(Tag).find('Compare').text
            ValueType    = XmlObj.getroot().find('QCFilters').find(Tag).find('CompareType').text

            FFunction = FilterFunction.FilterFunction(OperatorString=Operator,
                                                      CompareString=CompareValue,
                                                      CompareType=ValueType)

            InitLength  = len(DCs.DataContainers[ColumnTag].GetDataArray())
            FinalLength = None
            FilterArray = FFunction.Run(DataArray=DCs.DataContainers[ColumnTag].GetDataArray())

            if(not boDryRun):
                for Key in DCs.DataContainers.iterkeys():
                    DataArray = scipy.compress(FilterArray,
                                               DCs.DataContainers[Key].GetDataArray())
                    DCs.DataContainers[Key].ReplaceDataArray(DataArray)
                FinalLength = len(DCs.DataContainers[ColumnTag].GetDataArray())
            else:
                CounterDict = collections.defaultdict(int)
                for Entry in FilterArray:
                    CounterDict[Entry] += 1
                Difference  = CounterDict[False]
                FinalLength = InitLength-Difference

            self.SetFilterReportDictDict(ParentTag=ColumnTag,
                                         ChildTag=Tag,
                                         Value='Column Tag   = '+ColumnTag+'\n'+\
                                               'Filter Tag   = '+Tag+'\n'+\
                                               'Start Length = '+str(InitLength)+'\n'+\
                                               'Final Length = '+str(FinalLength)+'\n'+\
                                               'Difference   = '+str(InitLength-FinalLength))

        return DCs,\
               FilterTags

    def FilterAFCodedAlls(self,
                          XmlObj=lxml.etree._ElementTree,
                          DCs=DataContainer.DataContainers,
                          ColumnTag=str):

        FilterTags  = XmlObj.getroot().find('MtbGWAColumns').find(ColumnTag).find('Filters').text.split(',')

        self.InitFilterReportDictDict(ColumnTag)

        FilterArray = None
        for Tag in FilterTags:
            Operator     = XmlObj.getroot().find('QCFilters').find(Tag).find('Operator').text
            CompareValue = XmlObj.getroot().find('QCFilters').find(Tag).find('Compare').text
            ValueType    = XmlObj.getroot().find('QCFilters').find(Tag).find('CompareType').text

            FFunction = FilterFunction.FilterFunction(OperatorString=Operator,
                                                      CompareString=CompareValue,
                                                      CompareType=ValueType)

            InitLength  = len(DCs.DataContainers[ColumnTag].GetDataArray())
            FilterArray = FFunction.Run(DataArray=DCs.DataContainers[ColumnTag].GetDataArray())
            for Key in DCs.DataContainers.iterkeys():
                DataArray = scipy.compress(FilterArray,
                                           DCs.DataContainers[Key].GetDataArray())
                DCs.DataContainers[Key].ReplaceDataArray(DataArray)
            FinalLength = len(DCs.DataContainers[ColumnTag].GetDataArray())

            self.SetFilterReportDictDict(ParentTag=ColumnTag,
                                         ChildTag=Tag,
                                         Value='Column Tag   = '+ColumnTag+'\n'+\
                                               'Filter Tag   = '+Tag+'\n'+\
                                               'Start Length = '+str(InitLength)+'\n'+\
                                               'Final Length = '+str(FinalLength)+'\n'+\
                                               'Difference   = '+str(InitLength-FinalLength))

        return DCs,\
               FilterTags

    def FilterEMACs(self,
                    XmlObj=lxml.etree._ElementTree,
                    DCs=DataContainer.DataContainers,
                    ColumnTag=str):

        FilterTags  = XmlObj.getroot().find('MtbGWAColumns').find(ColumnTag).find('Filters').text.split(',')

        self.InitFilterReportDictDict(ColumnTag)

        FilterArray = None
        for Tag in FilterTags:
            Operator     = XmlObj.getroot().find('QCFilters').find(Tag).find('Operator').text
            CompareValue = XmlObj.getroot().find('QCFilters').find(Tag).find('Compare').text
            ValueType    = XmlObj.getroot().find('QCFilters').find(Tag).find('CompareType').text

            FFunction = FilterFunction.FilterFunction(OperatorString=Operator,
                                                      CompareString=CompareValue,
                                                      CompareType=ValueType)

            InitLength  = len(DCs.DataContainers[ColumnTag].GetDataArray())
            FilterArray = FFunction.Run(DataArray=DCs.DataContainers[ColumnTag].GetDataArray())
            for Key in DCs.DataContainers.iterkeys():
                DataArray = scipy.compress(FilterArray,
                                           DCs.DataContainers[Key].GetDataArray())
                DCs.DataContainers[Key].ReplaceDataArray(DataArray)
            FinalLength = len(DCs.DataContainers[ColumnTag].GetDataArray())

            self.SetFilterReportDictDict(ParentTag=ColumnTag,
                                         ChildTag=Tag,
                                         Value='Column Tag   = '+ColumnTag+'\n'+\
                                               'Filter Tag   = '+Tag+'\n'+\
                                               'Start Length = '+str(InitLength)+'\n'+\
                                               'Final Length = '+str(FinalLength)+'\n'+\
                                               'Difference   = '+str(InitLength-FinalLength))

        return DCs,\
               FilterTags

    def SetInitialArrayLength(self,
                              Value=int):
        self.InitialArrayLength = Value
        return

    def GetInitialArrayLength(self):
        return self.InitialArrayLength

    def SetNDeletedDuplicateSNPs(self,
                                 Value=int):
        self.NDeletedDuplicateSNPs = Value
        return

    def GetNDeletedDuplicateSNPs(self):
        return self.NDeletedDuplicateSNPs

    def RemoveDuplicateSNPs(self,
                            DCs=DataContainer.DataContainers):

        NRemoved   = 0
        DuplicateIndexDict = DCs.DataContainers['SNPID'].GetDuplicateIndexDict()
        for DCKey in DCs.DataContainers.iterkeys():
            NRemoved = DCs.DataContainers[DCKey].RemoveDuplicates(DuplicateIndexDict)

        self.SetNDeletedDuplicateSNPs(NRemoved)

        return DCs

    def SetMaxNDuplicateSNPs(self,
                             Value=int):
        self.MaxNDuplicateSNPs = Value
        return

    def GetMaxNDuplicateSNPs(self):
        return self.MaxNDuplicateSNPs

    def SetboDuplicateSNPWarning(self):
        if(self.GetMaxNDuplicateSNPs()>100):
            self.boDuplicateSNPWarning = True
        return

    def GetboDuplicateSNPWarning(self):
        return self.boDuplicateSNPWarning