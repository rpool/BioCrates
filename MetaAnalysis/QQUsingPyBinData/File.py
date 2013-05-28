import gzip
import re
import os
import scipy

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
        self.boUseLbzip2  = None
        self.boCompressed = None
        self.TmpDir       = None

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

    def ParseToDataContainers(self,
                              Delimiter=None):
#       Parse an input file into the DataContainers object
        DCs  = DataContainer.DataContainers()
        if(re.search('.npy',self.GetName())):
            Arrays = None
            if(self.GetboCompressed()):
                Arrays = scipy.load(self.GetDecomprName())
            else:
                Arrays = scipy.load(self.GetName())
            Header = Arrays[:,0].tolist()
            for i in xrange(len(Header)):
                Name                     = Header[i] # The names of the datacontainers are determined by the
                                                     # header column names.
                DCs.DataContainers[Name] = DataContainer.DataContainer()
                DCs.Names2Columns[Name]  = i
                DCs.Columns2Names[i]     = Name
                DCs.DataContainers[Name].SetDataArray(Arrays[i,1:])
                DCs.DataContainers[Name].SetDataName(Name)
            del Arrays
        else:
            Line = self.GetFileHandle().readline()
            if(self.GetboHeader()):
                Line  = re.sub('#','',Line)
                Names = Line.strip().split(Delimiter) # The file should be space or tab delimited!
                for i in range(len(Names)):
                    Name                     = Names[i] # The names of the datacontainers are determined by the
                                                        # header column names.
                    DCs.DataContainers[Name] = DataContainer.DataContainer()
                    DCs.Names2Columns[Name]  = i
                    DCs.Columns2Names[i]     = Name
                    DCs.DataContainers[Name].InitDataArray()
                    DCs.DataContainers[Name].SetDataName(Name)
            else:
                LSplit = Line.strip().split(Delimiter)
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
                LSplit = Line.strip().split(Delimiter)
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
        if(self.GetFileHandle()!=None):
            self.GetFileHandle().close()
        return

    def SetDecomprName(self):
#       Sets the name of the decompressed file, to be located in the temporary directory
        if(self.GetboCompressed()):
            if(re.search('.bz2',os.path.basename(self.GetName()))):
                self.DecomprName = re.sub('.bz2','',os.path.basename(self.GetName()))
            elif(re.search('.gz',os.path.basename(self.GetName()))):
                self.DecomprName = re.sub('.gz','',os.path.basename(self.GetName()))
            self.DecomprName = os.path.join(self.GetTmpDir(),self.DecomprName)
        return

    def GetDecomprName(self):
        return self.DecomprName

    def SetboCompressed(self):
#       Check if file is decompressed. For now only the .gz format!
        if(re.search('.gz',self.GetName())):
            self.boCompressed = True
        elif(re.search('.bz2',self.GetName())):
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

    def SetboUseLbzip2(self,
                       boUseLbzip2=bool):
        self.boUseLbzip2 = boUseLbzip2
        return

    def GetboUsePigz(self):
#       Pigz is the multi-threaded equivalent of gzip and gunzip => should go faster.
#       For now, there is no check on whether pigz is installed on the system.
        if(self.boUsePigz==None):
            self.boUsePigz = False # default value
        return self.boUsePigz

    def GetboUseLbzip2(self):
#       lbzip2 is the multi-threaded equivalent of lbzip2 => should go faster.
#       For now, there is no check on whether lbzip2 is installed on the system.
        if(self.boUseLbzip2==None):
            self.boUseLbzip2 = False # default value
        return self.boUseLbzip2

    def SetFileHandle(self,
                      Mode='r'): # 'r' to be safe :-)
        if(self.GetboCompressed()):
            if(self.GetboUsePigz()):
#               Decompress into the decompressed file located in the temporary directory.
                os.system('pigz -d -k -f -c '+self.GetName()+' > '+self.GetDecomprName())
                if(not re.search('.npy',self.GetDecomprName())):
                    self.FileHandle = open(self.GetDecomprName(),Mode)
            if(self.GetboUseLbzip2()):
#               Decompress into the decompressed file located in the temporary directory.
                os.system('lbzip2 -d -k -f -c '+self.GetName()+' > '+self.GetDecomprName())
                if(not re.search('.npy',self.GetDecomprName())):
                    self.FileHandle = open(self.GetDecomprName(),Mode)
            else:
                if(not re.search('.npy',self.GetName())):
                    self.FileHandle = gzip.open(self.GetName(),Mode)
        else:
            if(not re.search('.npy',self.GetName())):
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