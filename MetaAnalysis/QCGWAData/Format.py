import os
import re
import scipy
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
                    Path=str,
                    FilePreExtName=str):
        LogString = '  ++ Checking format of extra info files ...'
        print LogString
        Log.Write(LogString+'\n')

        FName  = os.path.join(Path,FilePreExtName)
        FmtLog = Logger.Logger(FileName=FName,
                               Extension='')
        LogString  = '## START TIMESTAMP\n'
        LogString += str(FmtLog.GetStartDate())+'\n'
        LogString += '## END TIMESTAMP'
        FmtLog.Write(LogString+'\n')
        LogString = Log.GetStartLogString()
        FmtLog.Write(LogString+'\n')

        self.CheckColumnFormat(DCsDict=DCsDict,
                               Log=Log,
                               FmtLog=FmtLog)
        self.CheckSNPIDs(DCsDict=DCsDict,
                         Log=Log,
                         FmtLog=FmtLog)
        self.CheckChrs(DCsDict=DCsDict,
                       Log=Log,
                       FmtLog=FmtLog)
        self.CheckStrandGenomes(DCsDict=DCsDict,
                                Log=Log,
                                FmtLog=FmtLog)
        self.CheckImputeds(DCsDict=DCsDict,
                           Log=Log,
                           FmtLog=FmtLog)
        self.CheckUsedForImps(DCsDict=DCsDict,
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

    def CheckUsedForImps(self,
                         DCsDict={},
                         Log=Logger,
                         FmtLog=Logger):
        for FName, DCs in DCsDict.iteritems():
            LogString = '    ++ Checking \"used_for_imp\" fields format for \"'+FName+'\" ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

            FilterArray    = (DCs.DataContainers['used_for_imp'].GetDataArray()=='1')
            TmpFilterArray = (DCs.DataContainers['used_for_imp'].GetDataArray()=='0')
            FilterArray    = (FilterArray | TmpFilterArray)
            LogString = '    ** Array \"used_for_imp\" complies with the ENGAGE analysis plan v3.0!'
            if(len(scipy.compress(FilterArray,DCs.DataContainers['imputed'].GetDataArray()))!=
               len(DCs.DataContainers['imputed'].GetDataArray())):
                LogString = '    ** Array \"used_for_imp\" does not comply with the ENGAGE analysis plan v3.0!'

            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')
            LogString = '    -- Done ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

        return

    def CheckImputeds(self,
                      DCsDict={},
                      Log=Logger,
                      FmtLog=Logger):
        for FName, DCs in DCsDict.iteritems():
            LogString = '    ++ Checking \"imputed\" fields format for \"'+FName+'\" ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

            FilterArray    = (DCs.DataContainers['imputed'].GetDataArray()=='1')
            TmpFilterArray = (DCs.DataContainers['imputed'].GetDataArray()=='0')
            FilterArray    = (FilterArray | TmpFilterArray)
            LogString = '    ** Array \"imputed\" complies with the ENGAGE analysis plan v3.0!'
            if(len(scipy.compress(FilterArray,DCs.DataContainers['imputed'].GetDataArray()))!=
               len(DCs.DataContainers['imputed'].GetDataArray())):
                LogString = '    ** Array \"imputed\" does not comply with the ENGAGE analysis plan v3.0!'

            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')
            LogString = '    -- Done ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

        return

    def CheckStrandGenomes(self,
                           DCsDict={},
                           Log=Logger,
                           FmtLog=Logger):
        for FName, DCs in DCsDict.iteritems():
            LogString = '    ++ Checking \"strand_genome\" fields format for \"'+FName+'\" ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

            FilterArray    = (DCs.DataContainers['strand_genome'].GetDataArray()=='+')
            TmpFilterArray = (DCs.DataContainers['strand_genome'].GetDataArray()=='-')
            FilterArray    = (FilterArray | TmpFilterArray)
            LogString = '    ** Array \"strand_genome\" complies with the ENGAGE analysis plan v3.0!'
            if(len(scipy.compress(FilterArray,DCs.DataContainers['strand_genome'].GetDataArray()))!=
               len(DCs.DataContainers['strand_genome'].GetDataArray())):
                LogString = '    ** Array \"strand_genome\" does not comply with the ENGAGE analysis plan v3.0!'

            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')
            LogString = '    -- Done ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

        return

    def CheckChrs(self,
                  DCsDict={},
                  Log=Logger,
                  FmtLog=Logger):
        for FName, DCs in DCsDict.iteritems():
            LogString = '    ++ Checking \"chr\" fields format for \"'+FName+'\" ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

            ChrRange = []
            for i in range(22):
                ChrRange.append(str(i+1))
            FilterArray = (DCs.DataContainers['chr'].GetDataArray()==ChrRange[0])
            for i  in range(1,len(ChrRange)):
                TmpFilterArray = (DCs.DataContainers['chr'].GetDataArray()==ChrRange[i])
                FilterArray    = (FilterArray | TmpFilterArray)
            LogString = '    ** Array \"chr\" complies with the ENGAGE analysis plan v3.0!'
            if(len(scipy.compress(FilterArray,DCs.DataContainers['chr'].GetDataArray()))!=
               len(DCs.DataContainers['chr'].GetDataArray())):
                LogString = '    ** Array \"chr\" does not comply with the ENGAGE analysis plan v3.0!'

            print LogString
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
            LogString = '    ++ Checking \"SNPID\" fields format for \"'+FName+'\" ...'
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
                LogString += '    ** Array \"SNPID\" complies with the ENGAGE analysis plan v3.0!'

            print LogString
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
