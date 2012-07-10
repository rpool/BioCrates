import os
import lxml.etree
import sys
import scipy
import re

import Logger
import File
import DataContainer

class HapMap:
    def __init__(self):
        self.SourceBuild                = None
        self.SourceRelease              = None
        self.SourceRefPanel             = None
        self.DestBuild                  = None
        self.DestRelease                = None
        self.DestRefPanel               = None
        self.DestFileName               = None
        self.DestFile                   = None
        self.DFile                      = None
        self.DestPath                   = None
        self.Delimiter                  = None
        self.Splitter                   = None
        self.LineArray                  = None
        self.ChrDict                    = {}
        self.PosDict                    = {}
        self.AllDict                    = {}
        self.StrandDict                 = {}
        self.MAFDict                    = {}
        self.NonHapMapIndexList         = None
        self.HapMapIndexList            = None
        self.StrandErrIndexList         = None
        self.AllErrIndexList            = None
        self.HapMapFalseImputedTrueList = None
        self.HapMapTrueImputedFalseList = None

        return

    def Clean(self):
        self.ResetNonHapMapIndexList()
        self.ResetHapMapIndexList()
        self.ResetStrandErrIndexList()
        self.ResetAllErrIndexList()
        self.ResetHapMapFalseImputedTrueList()
        self.ResetHapMapTrueImputedFalseList()
        return

    def SetNonHapMapIndexList(self,
                              List=[]):
        self.NonHapMapIndexList = List
        return

    def ResetNonHapMapIndexList(self):
        self.NonHapMapIndexList = None
        return

    def GetNonHapMapIndexList(self):
        return self.NonHapMapIndexList

    def SetHapMapIndexList(self,
                           List=[]):
        self.HapMapIndexList = List
        return

    def ResetHapMapIndexList(self):
        self.HapMapIndexList = None
        return

    def GetHapMapIndexList(self):
        return self.HapMapIndexList

    def SetStrandErrIndexList(self,
                              List=[]):
        self.StrandErrIndexList = List
        return

    def ResetStrandErrIndexList(self):
        self.StrandErrIndexList = None
        return

    def GetStrandErrIndexList(self):
        return self.StrandErrIndexList

    def SetAllErrIndexList(self,
                           List=[]):
        self.AllErrIndexList = List
        return

    def ResetAllErrIndexList(self):
        self.AllErrIndexList = None
        return

    def GetAllErrIndexList(self):
        return self.AllErrIndexList

    def SetHapMapFalseImputedTrueList(self,
                                      List=[]):
        self.HapMapFalseImputedTrueList = List
        return

    def ResetHapMapFalseImputedTrueList(self):
        self.HapMapFalseImputedTrueList = None
        return

    def GetHapMapFalseImputedTrueList(self):
        return self.HapMapFalseImputedTrueList

    def SetHapMapTrueImputedFalseList(self,
                                      List=[]):
        self.HapMapTrueImputedFalseList = List
        return

    def ResetHapMapTrueImputedFalseList(self):
        self.HapMapTrueImputedFalseList = None
        return

    def GetHapMapTrueImputedFalseList(self):
        return self.HapMapTrueImputedFalseList

    def Report(self,
               ReportName=str,
               SNPIDArray=scipy.array):

        ReportPath = os.getcwd()
        ReportPath = os.path.join(ReportPath,'HapMap')
        if(not os.path.isdir(ReportPath)):
            os.mkdir(ReportPath)
        ReportFile = os.path.join(ReportPath,ReportName)

        Ext    = ''
        Report = Logger.Logger(ReportFile,
                               Ext)
        LogString  = '## START TIMESTAMP\n'
        LogString += str(Report.GetStartDate())+'\n'
        LogString += '## END TIMESTAMP'
        Report.Write(LogString+'\n')
        LogString = Report.GetStartLogString()
        Report.Write(LogString+'\n')

        LogString = 'Belows follows the report of the HapMap mapping procedure from build_release_panel='+\
                    self.GetSourceBuild()+\
                    '_'+\
                    self.GetSourceRelease()+\
                    '_'+\
                    self.GetSourceRefPanel()+\
                    ' to build_release_panel='+\
                    self.GetDestBuild()+\
                    '_'+\
                    self.GetDestRelease()+\
                    '_'+\
                    self.GetDestRefPanel()+\
                    '.'
        Report.Write(LogString+'\n')
        Report.Write('\n')
        LogString  = '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'
        LogString += '! Each SNP reported in this file has NOT been mapped to the destination HapMap build_release_panel !\n'
        LogString += '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        Report.Write(LogString+'\n')
        Report.Write('\n')

        if(len(self.GetNonHapMapIndexList())>0):
            LogString = '** Found the following non-HapMap rsids (N='+\
                        str(len(self.GetNonHapMapIndexList()))+\
                        '):'
            Report.Write(LogString+'\n')
            for Index in self.GetNonHapMapIndexList():
                LogString = '   '+SNPIDArray[Index]
                Report.Write(LogString+'\n')
            Report.Write('\n')
            if(len(self.GetHapMapFalseImputedTrueList())>0):
                LogString = '**** Of these, the following rsids WERE imputed for the source HapMap Build/Release (N='+\
                            str(len(self.GetHapMapFalseImputedTrueList()))+\
                            '):'
                Report.Write(LogString+'\n')
                for Index in self.GetHapMapFalseImputedTrueList():
                    LogString = '     '+SNPIDArray[Index]
                    Report.Write(LogString+'\n')
                Report.Write('\n')
        if(len(self.GetStrandErrIndexList())>0):
            LogString = '** The following rsids have a strand_genome mismatch between the GWA data (-) and HapMap(+) (N='+\
                        str(len(self.GetStrandErrIndexList()))+\
                        '):'
            Report.Write(LogString+'\n')
            for Index in self.GetStrandErrIndexList():
                LogString = '   '+SNPIDArray[Index]
                Report.Write(LogString+'\n')
            Report.Write('\n')
        if(len(self.GetAllErrIndexList())>0):
            LogString = '** The following rsids have an allele naming mismatch between the GWA data and HapMap (N='+\
                        str(len(self.GetAllErrIndexList()))+\
                        '):'
            Report.Write(LogString+'\n')
            for Index in self.GetAllErrIndexList():
                LogString = '   '+SNPIDArray[Index]
                Report.Write(LogString+'\n')
            Report.Write('\n')
        if(len(self.GetHapMapIndexList())>0):
            if(len(self.GetHapMapTrueImputedFalseList())>0):
                LogString = 'The following rsids WERE NOT imputed although they are listed in the destination HapMap Build/Release (N='+\
                            str(len(self.GetHapMapTrueImputedFalseList()))+\
                            '):'
                Report.Write(LogString+'\n')
                for Index in self.GetAllErrIndexList():
                    LogString = '   '+SNPIDArray[Index]
                    Report.Write(LogString+'\n')
                Report.Write('\n')

        LogString = Report.GetEndLogString()
        Report.Write(LogString+'\n')
        Report.Close()
        return

    def Map(self,
            DCs=DataContainer.DataContainers):
        HapMapIndexList    = []
        NonHapMapIndexList = []
        StrandErrIndexList = []
        AllErrIndexList    = []
        HapMapMAFDataArray = []
        HapMapFalseImputedTrueList = []
        HapMapTrueImputedFalseList = []
        for i in range(len(DCs.DataContainers['SNPID'].GetDataArray())):
            SNPID = DCs.DataContainers['SNPID'].GetDataArray()[i]
            if(self.ChrDict.has_key(SNPID)):
                HapMapIndexList.append(i)
                HapMapMAFDataArray.append(self.MAFDict[SNPID])
                boDoMap = True
                if(DCs.DataContainers['strand_genome'].GetDataArray()[i]!=self.StrandDict[SNPID]):
                    StrandErrIndexList.append(i)
                    boDoMap = False
                if((not DCs.DataContainers['coded_all'].GetDataArray()[i] in self.AllDict[SNPID]) and
                   (not DCs.DataContainers['noncoded_all'].GetDataArray()[i] in self.AllDict[SNPID])):
                    AllErrIndexList.append(i)
                    boDoMap = False
                if(DCs.DataContainers['imputed'].GetDataArray()[i]=='0'):
                    HapMapTrueImputedFalseList.append(i)
                if(boDoMap):
                    DCs.DataContainers['position'].ReplaceDataArrayEntryAtIndex(Index=i,
                                                                                Value=self.PosDict[SNPID])
                    DCs.DataContainers['chr'].ReplaceDataArrayEntryAtIndex(Index=i,
                                                                           Value=self.ChrDict[SNPID])
            else:
                NonHapMapIndexList.append(i)
                if(DCs.DataContainers['imputed'].GetDataArray()[i]=='0'):
                    HapMapFalseImputedTrueList.append(i)
                HapMapMAFDataArray.append('NA')

        self.SetNonHapMapIndexList(List=NonHapMapIndexList)
        self.SetHapMapIndexList(List=HapMapIndexList)
        self.SetStrandErrIndexList(List=StrandErrIndexList)
        self.SetAllErrIndexList(List=AllErrIndexList)
        self.SetHapMapFalseImputedTrueList(List=HapMapFalseImputedTrueList)
        self.SetHapMapTrueImputedFalseList(List=HapMapTrueImputedFalseList)

        return DCs

    def ProcessLineArray(self):

        for i in range(1,len(self.GetLineArray())): # skip header
            Line   = self.GetLineArray()[i]

            NAll   = 2
            if(re.search('NAlleleNames=',Line)):
                NAll = Line.split('=')[-1]
                NAll = NAll.split('>')[0].strip()
                NAll = int(NAll)

            LSplit = Line.split(self.GetSplitter())
            SNPID  = LSplit[0]
            Chr    = LSplit[1]
            Pos    = LSplit[2]
            Alls   = LSplit[3:3+NAll]
            Strand = LSplit[3+NAll]
            MAF    = LSplit[3+NAll+1]

            self.ChrDict[SNPID]    = Chr
            self.PosDict[SNPID]    = Pos
            self.AllDict[SNPID]    = Alls
            self.StrandDict[SNPID] = Strand
            self.MAFDict[SNPID]    = MAF

        return


    def CheckIfFileExists(self,
                          Log=Logger,
                          HeadingSpaces=''):
        LogString  = HeadingSpaces
        LogString += '  ++ Checking if \"'+self.GetDestFile()+'\" exist(s) and are (is a) regular file(s) ...'
        print LogString
        Log.Write(LogString+'\n')

        if(os.path.isfile(self.GetDestFile()) or
           os.path.islink(self.GetDestFile())):
            LogString  = HeadingSpaces
            LogString += '    ** Path \"'+self.GetDestFile()+'\" is a regular file or a symlink to one ...'
            print LogString
            Log.Write(LogString+'\n')
        else:
            LogString  = HeadingSpaces
            LogString += '    ** Path \"'+self.GetDestFile()+'\" is not a regular file or a symlink to one...\n'
            LogString += '!!! EXITING !!!'
            print LogString
            Log.Write(LogString+'\n')
            sys.exit(1)

        LogString  = HeadingSpaces
        LogString += '  -- Done ...'
        print LogString
        Log.Write(LogString+'\n')

        return

    def SetDelimiter(self,
                     Value=str):
        self.Delimiter = Value
        return

    def GetDelimiter(self):
        return self.Delimiter

    def SetSplitter(self):

        if(self.GetDelimiter()=='WhiteSpace'):
            self.Splitter = None
        elif(self.GetDelimiter()==','):
            self.Splitter = ','

        return

    def GetSplitter(self):
        return self.Splitter

    def SetLineArray(self,
                     Array=scipy.array):
        self.LineArray = Array
        return

    def GetLineArray(self):
        return self.LineArray

    def ParseDestFile(self,
                      Log=Logger,
                      HeadingSpaces=''):

        LogString  = HeadingSpaces
        LogString += '  ++ Constructing data structure for \"'+self.GetDestFile()+'\" ...'
        print LogString
        Log.Write(LogString+'\n')

        self.DFile = File.File(Name=self.GetDestFile(),
                               boHeader=True)
        self.DFile.SetboUsePigz(boUsePigz=True)
        self.DFile.SetFileHandle(Mode='r')

        NLinesInFile,\
        NLinesInArray  = self.DFile.ParseToLineArray()
        self.SetLineArray(self.DFile.GetLineArray())

        LogString  = HeadingSpaces
        LogString += '    ** Removed '+str(NLinesInFile-NLinesInArray)+' duplicate lines!'
        print LogString
        Log.Write(LogString+'\n')

        self.DFile.Close()
        self.DFile.Cleanup()

        LogString  = HeadingSpaces
        LogString += '  -- Done ...'
        print LogString
        Log.Write(LogString+'\n')

        return

    def ProcessXml(self,
                   XmlObj=lxml.etree._ElementTree,
                   Tag=str,
                   Log=Logger):
        HapMapXml = XmlObj.getroot().find(Tag)
        self.SetSourceBuild(HapMapXml.find('SourceBuild').text)
        self.SetSourceRelease(HapMapXml.find('SourceRelease').text)
        self.SetSourceRefPanel(HapMapXml.find('SourceRefPanel').text)
        self.SetDestBuild(HapMapXml.find('DestBuild').text)
        self.SetDestRelease(HapMapXml.find('DestRelease').text)
        self.SetDestRefPanel(HapMapXml.find('DestRefPanel').text)
        self.SetDestFileName(HapMapXml.find('SummaryFile').text)
        self.SetDestPath(HapMapXml.find('SummaryPath').text)
        self.SetDestFile(os.path.join(self.GetDestPath(),self.GetDestFileName()))
        self.SetDelimiter(HapMapXml.find('Delimiter').text)
        self.SetSplitter()

        return

    def SetDestFileName(self,
                        Value=str):
        self.DestFileName = Value
        return

    def GetDestFileName(self):
        return self.DestFileName

    def SetDestFile(self,
                    Value=str):
        self.DestFile = Value
        return

    def GetDestFile(self):
        return self.DestFile

    def SetDestPath(self,
                    Value=str):
        self.DestPath = Value
        return

    def GetDestPath(self):
        return self.DestPath

    def SetSourceBuild(self,
                       Value=str):
        self.SourceBuild = Value
        return

    def GetSourceBuild(self):
        return self.SourceBuild

    def SetSourceRelease(self,
                         Value=str):
        self.SourceRelease = Value
        return

    def GetSourceRelease(self):
        return self.SourceRelease

    def SetSourceRefPanel(self,
                          Value=str):
        self.SourceRefPanel = Value
        return

    def GetSourceRefPanel(self):
        return self.SourceRefPanel

    def SetDestBuild(self,
                     Value=str):
        self.DestBuild = Value
        return

    def GetDestBuild(self):
        return self.DestBuild

    def SetDestRelease(self,
                       Value=str):
        self.DestRelease = Value
        return

    def GetDestRelease(self):
        return self.DestRelease

    def SetDestRefPanel(self,
                        Value=str):
        self.DestRefPanel = Value
        return

    def GetDestRefPanel(self):
        return self.DestRefPanel

