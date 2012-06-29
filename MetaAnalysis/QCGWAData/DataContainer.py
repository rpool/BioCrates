import re
import gzip
import scipy
import sys
import collections

import Logger

#===============================================================================
# This module contains the basic DataContainer and DataContainers classes.
# Their members and member modules should speak for themselves, if not a
# comment is provided.
#===============================================================================

class DataContainer:
    def __init__(self):
        self.DataName           = None
        self.DataArray          = None
        self.DuplicateDict      = None
        self.DuplicateIndexDict = None
        return

    def SetDataName(self,
                    Name=str):
        self.DataName = Name
        return

    def GetDataName(self):
        return self.DataName

    def InitDataArray(self):
        self.DataArray = []
        return

    def AppendToArray(self,
                      Entry=str):
        self.DataArray.append(Entry)
        return

    def CastDataArrayToScipy(self):
        self.DataArray = scipy.array(self.GetDataArray())

    def GetDataArray(self):
        return self.DataArray

    def RenameFieldsInDataArray(self,
                                Source=str,
                                Dest=str):
        DataArray      = self.GetDataArray()
        ExtractedArray = scipy.where(DataArray==Source)[0]
        for Index in ExtractedArray:
            DataArray[Index] = Dest
        self.DataArray = DataArray
        return

    def RenameColumnOfDataArray(self,
                                Source,
                                Dest):
        self.SetDataName(Dest)
        return

    def FindDuplicates(self):
        self.CounterDict   = collections.defaultdict(int)
        DuplicateIndexDict = collections.defaultdict(list)
        for i in range(len(self.GetDataArray())):
            Entry = self.GetDataArray()[i]
            self.CounterDict[Entry] += 1
            DuplicateIndexDict[Entry].append(i)
        self.DuplicateDict      = {}
        self.DuplicateIndexDict = {}
        for Key, Value in self.CounterDict.iteritems():
            if(Value>1):
                self.DuplicateDict[Key]      = Value
                self.DuplicateIndexDict[Key] = DuplicateIndexDict[Key]
        return self.DuplicateDict

    def GetDuplicateIndexDict(self):
        return self.DuplicateIndexDict

    def RemoveDuplicates(self,
                         DuplicateIndexDict={}):
        DataArray = self.GetDataArray()
        DelList   = []
        for Key in DuplicateIndexDict.iterkeys():
            for i in range(1,len(DuplicateIndexDict[Key])):
                Index = DuplicateIndexDict[Key][i]
                DelList.append(Index)
        DelList.sort()
        DelList.reverse()
        DataArray = scipy.delete(DataArray,tuple(DelList))

        self.DataArray = DataArray
        return

    def InitEntry2IndexDict(self):
        self.Entry2IndexDict = {}
        return

    def SetEntry2IndexDict(self):
        if(self.GetEntry2IndexDict()==None):
            self.InitEntry2IndexDict()
        for i in range(len(self.GetDataArray())):
            Entry = self.GetDataArray()[i]
            self.Entry2IndexDict[Entry] = i
        return

    def GetEntry2IndexDict(self):
        return self.Entry2IndexDict

class DataContainers:
    def __init__(self):
        self.DataContainers = {}
        self.Names2Columns  = {}
        self.Columns2Names  = {}
        return