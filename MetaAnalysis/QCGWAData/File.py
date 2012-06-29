import gzip
import re
import os
import lxml.etree
import scipy
import collections

import DataContainer

#===============================================================================
# This module contains the basic File class.
# Its members and member modules should speak for themselves, if not a
# comment is provided.
#===============================================================================

class File:
    def __init__(self,
                 Name=str,
                 boHeader=False):
        self.Name         = None
        self.DecomprName  = None
        self.boHeader     = None
        self.FileHandle   = None
        self.boUsePigz    = None
        self.boCompressed = None
        self.TmpDir       = None
        self.LineArray    = None

        self.SetName(Name=Name)
        self.SetboHeader(boHeader=boHeader)
        self.SetboCompressed()
        if(self.GetboCompressed()):
            self.SetTmpDir()
            self.SetDecomprName()
        return

    def SetTmpDir(self):
#       Sets the temporary directory
        self.TmpDir = os.path.join(os.getcwd(),'Tmp')
        if(not os.path.isdir(self.TmpDir)):
            os.mkdir(self.TmpDir)
        return

    def GetTmpDir(self):
        return self.TmpDir

    def FindDuplicates(self,
                       Array=scipy.array):
        CounterDict           = collections.defaultdict(int)
        TmpDuplicateIndexDict = collections.defaultdict(list)
        for i in range(len(Array)):
            Entry               = Array[i]
            CounterDict[Entry] += 1
            TmpDuplicateIndexDict[Entry].append(i)
        DuplicateDict      = {}
        DuplicateIndexDict = {}
        for Key, Value in CounterDict.iteritems():
            if(Value>1):
                DuplicateDict[Key]      = Value
                DuplicateIndexDict[Key] = TmpDuplicateIndexDict[Key]
        return DuplicateDict,\
               DuplicateIndexDict

    def RemoveDuplicates(self,
                         DuplicateIndexDict={},
                         DataArray=scipy.array):
        DelList   = []
        for Key in DuplicateIndexDict.iterkeys():
            for i in range(1,len(DuplicateIndexDict[Key])):
                Index = DuplicateIndexDict[Key][i]
                DelList.append(Index)
        DelList.sort()
        DelList.reverse()
        DataArray = scipy.delete(DataArray,tuple(DelList))

        return DataArray

    def ParseToLineArray(self):
        LineArray = []
        for Line in self.GetFileHandle():
            LJoin = Line.strip().split()
            LJoin = ' '.join(LJoin)
            LineArray.append(LJoin)
        LineArray      = scipy.array(LineArray)
        DuplicateDict,\
        DuplicateIndexDict = self.FindDuplicates(LineArray)
        self.LineArray     = self.RemoveDuplicates(DuplicateIndexDict,
                                                   LineArray)

        return len(LineArray),len(self.LineArray)

    def GetLineArray(self):
        return self.LineArray

    def LineArray2DataContainers(self):
        DCs  = DataContainer.DataContainers()
        if(self.GetboHeader()):
            Header = re.sub('#','',self.GetLineArray()[0])
            Names  = Header.strip().split()
            for i in range(len(Names)):
                Name                     = Names[i] # The names of the datacontainers are determined by the
                                                    # header column names.
                DCs.DataContainers[Name] = DataContainer.DataContainer()
                DCs.Names2Columns[Name]  = i
                DCs.Columns2Names[i]     = Name
                DCs.DataContainers[Name].InitDataArray()
                DCs.DataContainers[Name].SetDataName(Name)
        else:
            Line   = self.GetLineArray()[0]
            LSplit = Line.strip().split()
            for i in range(len(LSplit)):
                Name                     = str(i)
                DCs.DataContainers[Name] = DataContainer.DataContainer()
                DCs.Names2Columns[Name]  = i
                DCs.Columns2Names[i]     = Name
                DCs.DataContainers[Name].InitDataArray()
                DCs.DataContainers[Name].SetDataName(Name)
                Entry = LSplit[i]
                DCs.DataContainers[Name].AppendToArray(Entry)

        for i in range(1,len(self.GetLineArray())):
            Line = self.GetLineArray()[i]
            LSplit = Line.strip().split()
            for i in range(len(LSplit)):
                Name  = DCs.Columns2Names[i]
                Entry = LSplit[i]
                DCs.DataContainers[Name].AppendToArray(Entry)

        del self.LineArray
        self.LineArray = None

        for Key in DCs.DataContainers.iterkeys():
            DCs.DataContainers[Key].CastDataArrayToScipy() # Make scipy.arrays of the lists.

        return DCs

    def ParseAsXml(self):
        XmlObject = lxml.etree.parse(self.GetName())
        return XmlObject

    def ParseToDataContainers(self):
#       Parse an input file into the DataContainers object
        DCs  = DataContainer.DataContainers()
        Line = self.GetFileHandle().readline()
        if(self.GetboHeader()):
            Line  = re.sub('#','',Line)
            Names = Line.strip().split() # The file should be space or tab delimited!
            for i in range(len(Names)):
                Name                     = Names[i] # The names of the datacontainers are determined by the
                                                    # header column names.
                DCs.DataContainers[Name] = DataContainer.DataContainer()
                DCs.Names2Columns[Name]  = i
                DCs.Columns2Names[i]     = Name
                DCs.DataContainers[Name].InitDataArray()
                DCs.DataContainers[Name].SetDataName(Name)
        else:
            LSplit = Line.strip().split()
            for i in range(len(LSplit)):
                Name                     = str(i)
                DCs.DataContainers[Name] = DataContainer.DataContainer()
                DCs.Names2Columns[Name]  = i
                DCs.Columns2Names[i]     = Name
                DCs.DataContainers[Name].InitDataArray()
                DCs.DataContainers[Name].SetDataName(Name)
                Entry = LSplit[i]
                DCs.DataContainers[Name].AppendToArray(Entry)

        for Line in self.GetFileHandle():
            LSplit = Line.strip().split()
            for i in range(len(LSplit)):
                Name  = DCs.Columns2Names[i]
                Entry = LSplit[i]
                DCs.DataContainers[Name].AppendToArray(Entry)

        for Key in DCs.DataContainers.iterkeys():
            DCs.DataContainers[Key].CastDataArrayToScipy() # Make scipy.arrays of the lists.

        return DCs

    def Cleanup(self):
        if(self.GetboCompressed() and self.GetboUsePigz()):
#           Remove the temporary file
            os.remove(self.GetDecomprName())
        return

    def Close(self):
        self.GetFileHandle().close()
        return

    def SetDecomprName(self):
#       Sets the name of the decompressed file, to be located in the temporary directory
        if(self.GetboCompressed()):
            self.DecomprName = re.sub('.gz','',os.path.basename(self.GetName()))
            self.DecomprName = os.path.join(self.GetTmpDir(),self.DecomprName)
        return

    def GetDecomprName(self):
        return self.DecomprName

    def SetboCompressed(self):
#       Check if file is decompressed. For now only the .gz format!
        if(re.search('.gz',self.GetName())):
            self.boCompressed = True
        else:
            self.boCompressed = False
        return

    def GetboCompressed(self):
        return self.boCompressed

    def SetboUsePigz(self,
                     boUsePigz=bool):
        self.boUsePigz = boUsePigz
        return

    def GetboUsePigz(self):
#       Pigz is the multi-threaded equivalent of gzip and gunzip => should go faster.
#       For now, there is no check on whether pigz is installed on the system.
        if(self.boUsePigz==None):
            self.boUsePigz = False # default value
        return self.boUsePigz

    def SetFileHandle(self,
                      Mode='r'): # 'r' to be safe :-)
        if(self.GetboCompressed()):
            if(self.GetboUsePigz()):
#               Decompress into the decompressed file located in the temporary directory.
                os.system('pigz -d -k -f -c '+self.GetName()+' > '+self.GetDecomprName())
                self.FileHandle = open(self.GetDecomprName(),Mode)
            else:
                self.FileHandle = gzip.open(self.GetName(),Mode)
        else:
            self.FileHandle = open(self.GetName(),Mode)
        return

    def GetFileHandle(self):
        return self.FileHandle

    def SetboHeader(self,
                    boHeader=bool):
#       Sets the boHeader flag: True -> file has a header; False -> file has no header
        self.boHeader = boHeader
        return

    def GetboHeader(self):
        return self.boHeader

    def SetName(self,
                Name=str):
        self.Name = Name
        return

    def GetName(self):
        return self.Name