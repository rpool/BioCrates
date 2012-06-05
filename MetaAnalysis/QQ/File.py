import gzip
import re
import os

import DataContainer

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

        self.SetName(Name=Name)
        self.SetboHeader(boHeader=boHeader)
        self.SetboCompressed()
        self.SetDecomprName()
        return

    def ParseToDataContainers(self):
        DCs  = DataContainer.DataContainers()
        Line = self.GetFileHandle().readline()
        if(self.GetboHeader()):
            Line  = re.sub('#','',Line)
            Names = Line.strip().split()
            for i in range(len(Names)):
                Name                     = Names[i]
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
            DCs.DataContainers[Key].CastDataArrayToScipy()

        return DCs

    def Cleanup(self):
        if(self.GetboCompressed() and self.GetboUsePigz()):
            os.remove(self.GetDecomprName())
        return

    def Close(self):
        self.GetFileHandle().close()
        return

    def SetDecomprName(self):
        if(self.GetboCompressed()):
            self.DecomprName = re.sub('.gz','',self.GetName())
        return

    def GetDecomprName(self):
        return self.DecomprName

    def SetboCompressed(self):
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
        if(self.boUsePigz==None):
            self.boUsePigz = False # default value
        return self.boUsePigz

    def SetFileHandle(self,
                      Mode='r'): # 'r' to be safe :-)
        if(self.GetboCompressed()):
            if(self.GetboUsePigz()):
                os.system('pigz -d -k -f '+self.GetName())
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