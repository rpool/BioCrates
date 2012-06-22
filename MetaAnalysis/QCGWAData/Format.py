import os
import re
import lxml.etree

import File
import Logger
import DataContainer

class Format:
    def __init__(self):
        self.ExtraInfoFiles = None
        self.ColumnFormat   = None
        self.Delimiter      = None
        self.Split          = None
        return

    def SetExtraInfoFiles(self):
        self.ExtraInfoFiles = []
        return

    def GetExtraInfoFiles(self):
        if(self.ExtraInfoFiles==None):
            self.ExtraInfoFiles = []
        return self.ExtraInfoFiles

    def AppendFilesToExtraInfoFiles(self,
                                    XmlObj=lxml.etree._ElementTree,
                                    Log=Logger):
        for F in XmlObj.getroot().find('ExtraInfoFiles'):
            if(eval(F.find('boUse').text)):
                FName = os.path.join(F.find('Path').text,F.find('Name').text)

                LogString = '  ++ Constructing data structure for \"'+FName+'\" ...'
                print LogString
                Log.Write(LogString+'\n')

                FFile = File.File(Name=FName,
                                  boHeader=True)
                FFile.SetboUsePigz(boUsePigz=True)
                FFile.SetFileHandle(Mode='r')
                self.GetExtraInfoFiles().append(FFile)

                LogString = '  -- Done ...'
                print LogString
                Log.Write(LogString+'\n')
        return

    def ParseExtraInfoFiles(self,
                            Log=Logger):
        DCsDict = {}
        for EIFile in self.GetExtraInfoFiles():
            LogString = '  ++ Parsing \"'+EIFile.GetName()+'\" ...'
            print LogString
            Log.Write(LogString+'\n')

            DCsDict[EIFile.GetName()] = EIFile.ParseToDataContainers()
            EIFile.Close()
            EIFile.Cleanup()

            LogString = '  -- Done ...'
            print LogString
            Log.Write(LogString+'\n')
        return DCsDict

    def SetSplitFunction(self,
                         Log=Logger):
        FuncName = None
        if(self.GetDelimiter()=='WhiteSpace'):
            self.Splitter = None
            FuncName   = 'str.split(None)'
        elif(self.GetDelimiter()==','):
            self.Splitter = ','
            FuncName   = 'str.split(\',\')'
        LogString = '  ++ Set str.split() function according to \"'+FuncName+'\" ...'
        print LogString
        Log.Write(LogString+'\n')
        LogString = '  -- Done ...'
        print LogString
        Log.Write(LogString+'\n')
        return

    def SetDelimiter(self,
                     XmlObj=lxml.etree._ElementTree,
                     Log=Logger):
        self.Delimiter = XmlObj.getroot().find('Format').find('Delimiter').text
        LogString = '  ++ Set delimiter to '+self.Delimiter
        print LogString
        Log.Write(LogString+'\n')
        LogString = '  -- Done ...'
        print LogString
        Log.Write(LogString+'\n')
        return

    def GetDelimiter(self):
        if(self.Delimiter==None):
            self.Delimiter = 'WhiteSpace' # Set to default
        return self.Delimiter

    def SetColumnFormat(self,
                        XmlObj=lxml.etree._ElementTree,
                        Log=Logger):
        self.ColumnFormat = XmlObj.getroot().find('Format').find('ColumnNames').text.split(',')
        LogString = '  ++ Set column format list to ['+(',').join(self.ColumnFormat)+'] ...'
        print LogString
        Log.Write(LogString+'\n')
        LogString = '  -- Done ...'
        print LogString
        Log.Write(LogString+'\n')
        return

    def GetColumnFormat(self):
        return self.ColumnFormat

    def CheckFormat(self,
                    DCsDict={},
                    Log=Logger,
                    Path=str):
        LogString = '  ++ Checking format of extra info files ...'
        print LogString
        Log.Write(LogString+'\n')

        FName  = os.path.join(Path,'CheckFormatExtraInfoFiles')
        FmtLog = Logger.Logger(FileName=FName,
                               Extension='')
        LogString  = '## START TIMESTAMP\n'
        LogString += str(FmtLog.GetStartDate())+'\n'
        LogString += '## END TIMESTAMP'
        FmtLog.Write(LogString+'\n')

        self.CheckColumnFormat(DCsDict=DCsDict,
                               Log=Log,
                               FmtLog=FmtLog)
        self.CheckSNPIDs(DCsDict=DCsDict,
                         Log=Log,
                         FmtLog=FmtLog)
        self.Checkchrs(DCsDict=DCsDict,
                       Log=Log,
                       FmtLog=FmtLog)

        LogString = '\n**** Done :-)'
        FmtLog.Write(LogString+'\n')
        LogString = FmtLog.GetEndLogString()
        FmtLog.Write(LogString+'\n')
        FmtLog.Close()

        LogString = '  -- Done ...'
        print LogString
        Log.Write(LogString+'\n')

        return

    def Checkchrs(self,
                  DCsDict={},
                  Log=Logger,
                  FmtLog=Logger):
        for FName, DCs in DCsDict.iteritems():
            LogString = '    ++ Checking SNPID fields format for \"'+FName+'\" ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')
            LogString  = ''
            boComplies = True
            for Entry in DCs.DataContainers['SNPID'].GetDataArray():
                if(Entry[0:2]!='rs'):
                    LogString += '    ** SNPID does not comply with the ENGAGE analysis plan v3.0!\n'
                    LogString += '    ** SNPID '+Entry+' does not start with \"rs\"! '
                    boComplies = False
                    break
                if(not re.search('[0-9]',Entry)):
                    LogString += '    ** SNPID does not comply with the ENGAGE analysis plan v3.0!\n'
                    LogString += '    ** SNPID '+Entry+' does not have a number! '
                    boComplies = False
                    break
            if(boComplies):
                LogString += '    ** All SNPIDs comply the ENGAGE analysis plan v3.0!'
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')
            LogString = '    -- Done ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

        return

    def CheckSNPIDs(self,
                    DCsDict={},
                    Log=Logger,
                    FmtLog=Logger):
        for FName, DCs in DCsDict.iteritems():
            LogString = '    ++ Checking SNPID fields format for \"'+FName+'\" ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')
            LogString  = ''
            boComplies = True
            for Entry in DCs.DataContainers['SNPID'].GetDataArray():
                if(Entry[0:2]!='rs'):
                    LogString += '    ** SNPID does not comply with the ENGAGE analysis plan v3.0!\n'
                    LogString += '    ** SNPID '+Entry+' does not start with \"rs\"! '
                    boComplies = False
                    break
                if(not re.search('[0-9]',Entry)):
                    LogString += '    ** SNPID does not comply with the ENGAGE analysis plan v3.0!\n'
                    LogString += '    ** SNPID '+Entry+' does not have a number! '
                    boComplies = False
                    break
            if(boComplies):
                LogString += '    ** All SNPIDs comply the ENGAGE analysis plan v3.0!'
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')
            LogString = '    -- Done ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

        return

    def CheckColumnFormat(self,
                          DCsDict={},
                          Log=Logger,
                          FmtLog=Logger):
        for FName, DCs in DCsDict.iteritems():
            LogString = '    ++ Checking column format for \"'+FName+'\" ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')
            for Key in DCs.DataContainers.iterkeys():
                LogString = '    ** Column name \"'+Key+'\" '
                if((Key in self.GetColumnFormat()) and
                   (DCs.DataContainers[Key].GetDataName() in self.GetColumnFormat())):
                    LogString += 'complies with the ENGAGE analysis plan v3.0!'
                else:
                    LogString += 'does not comply with the ENGAGE analysis plan v3.0!'
                print LogString
                Log.Write(LogString+'\n')
                FmtLog.Write(LogString+'\n')
            LogString = '    -- Done ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

        return
