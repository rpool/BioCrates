# -*- coding: utf-8 -*-
"""
This module contains the basic DataContainer and DataContainers classes.

.. moduleauthor:: René Pool <r.pool@vu.nl>

"""
import re
import scipy
import collections
import os

import Logger

#===============================================================================
# This module contains the basic DataContainer and DataContainers classes.
# Their members and member modules should speak for themselves, if not a
# comment is provided.
#===============================================================================

class DataContainer:
    """
    Class **DataContainer**
    """
    def __init__(self):
        self.DataName           = None
        self.DataArray          = None
        self.DuplicateDict      = None
        self.DuplicateIndexDict = None
        self.MaxNDuplicates     = None
        self.Entry2IndexDict    = None
        return

    def GetMaxNDuplicates(self):
        """
        Returns the maximum number of duplicate lines found in an input file.
        """
        return self.MaxNDuplicates

    def SetDataName(self,
                    Name=str):
        """
        Sets the DataContainer name.
        """
        self.DataName = Name
        return

    def GetDataName(self):
        """
        Returns the DataContainer name.
        """
        return self.DataName

    def InitDataArray(self):
        """
        Initializes the DataArray
        """
        self.DataArray = []
        return

    def ReplaceDataArray(self,
                         DataArray=scipy.array):
        """
        Replaces the DataArray by the input array.

        :param DataArray: Array that replaces the old value
        :type DataArray: scipy.array
        :returns: nothing
        """
        self.DataArray = DataArray
        return

    def ReplaceDataArrayEntryAtIndex(self,
                                     Index=int,
                                     Value=str):
        """
        Replaces element *Index* of DataArray by the input *Value*.

        :param Index: Index of the element that will be replaced
        :type Index: int
        :param Value: New value of *DataArray[Index]*
        :type Value: str
        :returns: nothing
        """
        self.DataArray[Index] = Value
        return

    def AppendToArray(self,
                      Entry=str):
        """
        Appends an element to the DataArray

        :param Entry: Value of the appended element
        :type Entry: str
        :returns: nothing
        """
        self.DataArray.append(Entry)
        return

    def CastDataArrayToScipy(self):
        """
        Type casts the DataArray to type `scipy.array`.
        """
        self.DataArray = scipy.array(self.GetDataArray())

    def GetDataArray(self):
        """
        Returns the DataArray.
        """
        return self.DataArray

    def RenameFieldsInDataArray(self,
                                Source=str,
                                Dest=str):
        """
        Renames DataArray elements that are equal to *Source* to *Dest*

        :param Source: Source value
        :type Source: str
        :param Dest: Replacement value
        :type Dest: str
        :returns: nothing
        """
        DataArray      = self.GetDataArray()
        ExtractedArray = scipy.where(DataArray==Source)[0]
        for Index in ExtractedArray:
            DataArray[Index] = Dest
        self.DataArray = DataArray
        return

    def RenameColumnOfDataArray(self,
                                Source=str,
                                Dest=str):
        """
        Renames the DataContainer name from *Source* to *Dest*.

        :param Source: Source value
        :type Source: str
        :param Dest: Replacement value
        :type Dest: str
        :returns: nothing
        """
        self.SetDataName(Dest)
        return

    def FindDuplicates(self):
        """
        Finds duplicate entries in DataArray.

        :returns: A dictionary of duplicates
        """
        self.CounterDict   = collections.defaultdict(int)
        DuplicateIndexDict = collections.defaultdict(list)
        for i in range(len(self.GetDataArray())):
            Entry = self.GetDataArray()[i]
            self.CounterDict[Entry] += 1
            DuplicateIndexDict[Entry].append(i)
        self.DuplicateDict      = {}
        self.DuplicateIndexDict = {}
        self.MaxNDuplicates     = 0
        for Key, Value in self.CounterDict.iteritems():
            if(Value>1):
                self.DuplicateDict[Key]      = Value
                self.DuplicateIndexDict[Key] = DuplicateIndexDict[Key]
                self.MaxNDuplicates          = max(self.MaxNDuplicates,Value-1)
        return self.DuplicateDict

    def GetDuplicateIndexDict(self):
        """
        Returns the duplicate index list.
        """
        return self.DuplicateIndexDict

    def RemoveDuplicates(self,
                         DuplicateIndexDict={}):
        """
        Removes the duplicates listed in *DuplicateIndexDict*.

        :param DuplicateIndexDict: The dictionary of the indices of the duplicates
        :type DuplicateIndexDict: dict
        :returns: The number of removed lines (*int*)
        """
        DataArray = self.GetDataArray()
        NRemoved  = len(DataArray)
        DelList   = []
        for Key in DuplicateIndexDict.iterkeys():
            for i in range(1,len(DuplicateIndexDict[Key])):
                Index = DuplicateIndexDict[Key][i]
                DelList.append(Index)
        DelList.sort()
        DelList.reverse()
        DataArray = scipy.delete(DataArray,tuple(DelList))

        self.DataArray = DataArray
        NRemoved      -= len(DataArray)
        return NRemoved

    def InitEntry2IndexDict(self):
        """
        Initializes Entry2IndexDict.
        """
        self.Entry2IndexDict = {}
        return

    def SetEntry2IndexDict(self):
        """
        Sets Entry2IndexDict.
        """
        if(self.GetEntry2IndexDict()==None):
            self.InitEntry2IndexDict()
        for i in range(len(self.GetDataArray())):
            Entry = self.GetDataArray()[i]
            self.Entry2IndexDict[Entry] = i
        return

    def GetEntry2IndexDict(self):
        """
        Returns Entry2IndexDict.
        """
        return self.Entry2IndexDict

class DataContainers:
    def __init__(self):
        self.DataContainers = {}
        self.Names2Columns  = {}
        self.Columns2Names  = {}
        return

    def WriteBioCratesGWAOutput(self,
                                FileName=str,
                                OutPath=str,
                                HeaderList=[],
                                Header2ColumnDict={}):

        DelList = []
        for i in range(len(HeaderList)):
            Entry    = HeaderList[i]
            ColumnId = Header2ColumnDict[Entry]
            if(not self.DataContainers.has_key(ColumnId)):
                DelList.append(i)
        DelList.sort()
        DelList.reverse()
        for Index in DelList:
            del HeaderList[Index]


        if(not os.path.isdir(OutPath)):
            os.mkdir(OutPath)
        FilePath = os.path.join(OutPath,FileName)

        FH = open(FilePath,'w')

        ColumnWidthList = []
        MaxWidth        = 0
        for Entry in HeaderList:
            MaxWidth = max(len(Entry),MaxWidth)

        for Entry in HeaderList:
            ColumnWidthList.append(MaxWidth+6)

        for i in range(len(HeaderList)):
            Entry        = HeaderList[i]
            FormatString = '{0:>'+str(ColumnWidthList[i])+'}'
            FH.write(FormatString.format(Entry))
        FH.write('\n')

        ArrayLenth = len(self.DataContainers['SNPID'].GetDataArray())
        for i in range(ArrayLenth):
            for j in range(len(HeaderList)):
                Entry       = HeaderList[j]
                ColumnId    = Header2ColumnDict[Entry]

                String = str(self.DataContainers[ColumnId].GetDataArray()[i])
                if(String=='NA'):
                    FormatString = '{0:>'+str(ColumnWidthList[j])+'}'
                    String       = FormatString.format(String)
                elif(ColumnId=='beta'):
                    FormatString = '{0:>'+str(ColumnWidthList[j])+'.5e}'
                    String       = FormatString.format(float(String))
                elif(ColumnId=='SE'):
                    FormatString = '{0:>'+str(ColumnWidthList[j])+'.5e}'
                    String       = FormatString.format(float(String))
                elif(ColumnId=='pval'):
                    FormatString = '{0:>'+str(ColumnWidthList[j])+'.5e}'
                    String       = FormatString.format(float(String))
                elif(ColumnId=='PValWald'):
                    FormatString = '{0:>'+str(ColumnWidthList[j])+'.5e}'
                    String       = FormatString.format(float(String))
                elif(ColumnId=='AF_coded_all'):
                    FormatString = '{0:>'+str(ColumnWidthList[j])+'.5e}'
                    String       = FormatString.format(float(String))
                elif(ColumnId=='HWE_pval'):
                    FormatString = '{0:>'+str(ColumnWidthList[j])+'.5e}'
                    String       = FormatString.format(float(String))
                elif(ColumnId=='n_total'):
                    FormatString = '{0:>'+str(ColumnWidthList[j])+'d}'
                    String       = FormatString.format(int(round(float(String),0)))
                elif(ColumnId=='oevar_imp'):
                    FormatString = '{0:>'+str(ColumnWidthList[j])+'.5e}'
                    String       = FormatString.format(float(String))
                FormatString = '{0:>'+str(ColumnWidthList[j])+'}'
                FH.write(FormatString.format(String))
            FH.write('\n')
        FH.close()

        Cwd = os.getcwd()
        os.chdir(OutPath)
        os.system('pigz -f '+FileName)
        os.chdir(Cwd)

        return