import os
import re
import scipy
import lxml.etree
import sys

import File
import Logger
import DataContainer

class Format:
    def __init__(self):
        self.ExtraInfoFiles      = None
        self.ColumnFormat        = None
        self.Delimiter           = None
        self.Split               = None
        self.GWADataFileName     = None
        self.boColumnFormatOK    = None
        self.boFieldFormatOKDict = None
        return

    def InitboFieldFormatOKDict(self):
        self.boFieldFormatOKDict = {}
        return

    def SetboFieldFormatOKDict(self,
                               ColumnID=str,
                               Flag=bool):
        if(self.boFieldFormatOKDict==None):
            self.InitboFieldFormatOKDict()
        self.boFieldFormatOKDict[ColumnID] = Flag
        return

    def GetboFieldFormatOKDict(self):
        return self.boFieldFormatOKDict

    def SetboColumnFormatOK(self,
                            Flag=bool):
        self.boColumnFormatOK = Flag
        return

    def GetboColumnFormatOK(self):
        return self.boColumnFormatOK


    def SetGWADataFileName(self,
                       Name=str):
        self.GWADataFileName = Name
        return

    def GetGWADataFileName(self):
        return self.GWADataFileName

    def SetExtraInfoFiles(self):
        self.ExtraInfoFiles = []
        return

    def GetExtraInfoFiles(self):
        if(self.ExtraInfoFiles==None):
            self.ExtraInfoFiles = []
        return self.ExtraInfoFiles

    def AppendFilesToExtraInfoFiles(self,
                                    XmlObj=lxml.etree._ElementTree,
                                    Log=Logger,
                                    HeadingSpaces=''):
        for F in XmlObj.getroot().find('ExtraInfoFiles'):
            if(eval(F.find('boUse').text)):
                FName = os.path.join(F.find('Path').text,F.find('Name').text)

                LogString  = HeadingSpaces
                LogString += '  ++ Constructing data structure for \"'+FName+'\" ...'
                print LogString
                Log.Write(LogString+'\n')

                FFile = File.File(Name=FName,
                                  boHeader=True)
                FFile.SetboUsePigz(boUsePigz=True)
                FFile.SetFileHandle(Mode='r')
                self.GetExtraInfoFiles().append(FFile)

                LogString  = HeadingSpaces
                LogString += '  -- Done ...'
                print LogString
                Log.Write(LogString+'\n')
        return

    def ParseGWADataFile(self,
                         Log=Logger,
                         HeadingSpaces='',
                         boRemoveDuplicateLines=False):
        DCsDict = {}

        LogString  = HeadingSpaces
        LogString += '  ++ Parsing \"'+self.GetGWADataFileName()+'\" ...'
        print LogString
        Log.Write(LogString+'\n')

        LogString  = HeadingSpaces
        LogString += '    ++ Constructing data structure for \"'+self.GetGWADataFileName()+'\" ...'
        print LogString
        Log.Write(LogString+'\n')

        FFile = File.File(Name=self.GetGWADataFileName(),
                          boHeader=True)
        FFile.SetboUsePigz(boUsePigz=True)
        FFile.SetFileHandle(Mode='r')

        LogString  = HeadingSpaces
        LogString += '    -- Done ...'
        print LogString
        Log.Write(LogString+'\n')

        if(boRemoveDuplicateLines):
            NLinesInFile,\
            NLinesInArray  = FFile.ParseToLineArray()

            LogString  = HeadingSpaces
            LogString += '  ** Removed '+str(NLinesInFile-NLinesInArray)+' duplicate lines!'
            print LogString
            Log.Write(LogString+'\n')

            DCsDict['GWADataFile'] = FFile.LineArray2DataContainers()
        else:
            DCsDict['GWADataFile'] = FFile.ParseToDataContainers()

        FFile.Close()
        FFile.Cleanup()

        LogString  = HeadingSpaces
        LogString += '  -- Done ...'
        print LogString
        Log.Write(LogString+'\n')

        return DCsDict

    def ParseExtraInfoFiles(self,
                            Log=Logger):
        DCsDict = {}
        for EIFile in self.GetExtraInfoFiles():
            LogString = '    ++ Parsing \"'+EIFile.GetName()+'\" ...'
            print LogString
            Log.Write(LogString+'\n')

            LogString = '      ++ Constructing data structure for \"'+EIFile.GetName()+'\" ...'
            print LogString
            Log.Write(LogString+'\n')

            EFile = File.File(Name=EIFile.GetName(),
                              boHeader=True)
            EFile.SetboUsePigz(boUsePigz=True)
            EFile.SetFileHandle(Mode='r')

            LogString = '      -- Done ...'
            print LogString
            Log.Write(LogString+'\n')

            DCsDict[EIFile.GetName()] = EIFile.ParseToDataContainers()
            EIFile.Close()
            EIFile.Cleanup()

            LogString = '    -- Done ...'
            print LogString
            Log.Write(LogString+'\n')
        return DCsDict

    def SetSplitFunction(self,
                         Log=Logger,
                         HeadingSpaces=''):
        FuncName = None
        if(self.GetDelimiter()=='WhiteSpace'):
            self.Splitter = None
            FuncName   = 'str.split(None)'
        elif(self.GetDelimiter()==','):
            self.Splitter = ','
            FuncName   = 'str.split(\',\')'
        LogString  = HeadingSpaces
        LogString += '  ++ Set str.split() function according to \"'+FuncName+'\" ...'
        print LogString
        Log.Write(LogString+'\n')
        LogString  = HeadingSpaces
        LogString += '  -- Done ...'
        print LogString
        Log.Write(LogString+'\n')
        return

    def SetDelimiter(self,
                     XmlObj=lxml.etree._ElementTree,
                     Log=Logger,
                     HeadingSpaces=''):
        self.Delimiter = XmlObj.getroot().find('Format').find('Delimiter').text
        if(self.Delimiter=='WhiteSpace'):
            self.Delimiter = None
        LogString  = HeadingSpaces
        LogString += '  ++ Set delimiter to \"'+str(self.Delimiter)+'\" ...'
        print LogString
        Log.Write(LogString+'\n')
        LogString  = HeadingSpaces
        LogString += '  -- Done ...'
        print LogString
        Log.Write(LogString+'\n')
        return

    def GetDelimiter(self):
        if(self.Delimiter==None):
            self.Delimiter = 'WhiteSpace' # Set to default
        return self.Delimiter

    def SetColumnFormat(self,
                        XmlObj=lxml.etree._ElementTree,
                        Log=Logger,
                        HeadingSpaces=''):
        self.ColumnFormat = XmlObj.getroot().find('Format').find('ColumnNames').text.split(',')
        LogString  = HeadingSpaces
        LogString += '  ++ Set column format list to ['+(',').join(self.ColumnFormat)+'] ...'
        print LogString
        Log.Write(LogString+'\n')
        LogString  = HeadingSpaces
        LogString += '  -- Done ...'
        print LogString
        Log.Write(LogString+'\n')
        return

    def GetColumnFormat(self):
        return self.ColumnFormat

    def CheckFormat(self,
                    DCsDict={},
                    Log=Logger,
                    HeadingSpaces='',
                    Path=str,
                    FilePreExtName=str,
                    FileType=str,
                    XmlObj=lxml.etree._ElementTree,
                    Tag=str):
        LogString  = HeadingSpaces
        LogString += '  ++ Checking format of '+FileType+' ...'
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
                               FmtLog=FmtLog,
                               HeadingSpaces='    ')
        CondList = XmlObj.getroot().find(Tag).find('SNPID').find('MandatoryFieldEntries').text
        CondList = CondList.split(',')
        self.CheckSNPIDs(DCsDict=DCsDict,
                         Log=Log,
                         FmtLog=FmtLog,
                         HeadingSpaces='    ',
                         ConditionList=CondList)
        CondList = XmlObj.getroot().find(Tag).find('chr').find('MandatoryFieldEntries').text
        CondList = CondList.split(',')
        self.CheckChrs(DCsDict=DCsDict,
                       Log=Log,
                       FmtLog=FmtLog,
                       HeadingSpaces='    ',
                       ConditionList=CondList)
        CondList = XmlObj.getroot().find(Tag).find('strand_genome').find('MandatoryFieldEntries').text
        CondList = CondList.split(',')
        self.CheckStrandGenomes(DCsDict=DCsDict,
                                Log=Log,
                                FmtLog=FmtLog,
                                HeadingSpaces='    ',
                                ConditionList=CondList)
        CondList = XmlObj.getroot().find(Tag).find('imputed').find('MandatoryFieldEntries').text
        CondList = CondList.split(',')
        self.CheckImputeds(DCsDict=DCsDict,
                           Log=Log,
                           FmtLog=FmtLog,
                           HeadingSpaces='    ',
                           ConditionList=CondList)
        CondList = XmlObj.getroot().find(Tag).find('used_for_imp').find('MandatoryFieldEntries').text
        CondList = CondList.split(',')
        self.CheckUsedForImps(DCsDict=DCsDict,
                              Log=Log,
                              FmtLog=FmtLog,
                              HeadingSpaces='    ',
                              ConditionList=CondList)

        if(Tag=='MtbGWAColumns'):
            CondList = XmlObj.getroot().find(Tag).find('position').find('MandatoryFieldEntries').text
            CondList = CondList.split(',')
            self.CheckPositions(DCsDict=DCsDict,
                                Log=Log,
                                FmtLog=FmtLog,
                                HeadingSpaces='    ',
                                ConditionList=CondList)
            CondList = XmlObj.getroot().find(Tag).find('coded_all').find('MandatoryFieldEntries').text
            CondList = CondList.split(',')
            self.CheckCodedAlls(DCsDict=DCsDict,
                                Log=Log,
                                FmtLog=FmtLog,
                                HeadingSpaces='    ',
                                ConditionList=CondList)
            CondList = XmlObj.getroot().find(Tag).find('noncoded_all').find('MandatoryFieldEntries').text
            CondList = CondList.split(',')
            self.CheckNonCodedAlls(DCsDict=DCsDict,
                                   Log=Log,
                                   FmtLog=FmtLog,
                                   HeadingSpaces='    ',
                                   ConditionList=CondList)
            CondList = XmlObj.getroot().find(Tag).find('beta').find('MandatoryFieldEntries').text
            CondList = CondList.split(',')
            self.CheckBetas(DCsDict=DCsDict,
                            Log=Log,
                            FmtLog=FmtLog,
                            HeadingSpaces='    ',
                            ConditionList=CondList)
            CondList = XmlObj.getroot().find(Tag).find('SE').find('MandatoryFieldEntries').text
            CondList = CondList.split(',')
            self.CheckSEs(DCsDict=DCsDict,
                          Log=Log,
                          FmtLog=FmtLog,
                          HeadingSpaces='    ',
                          ConditionList=CondList)
            CondList = XmlObj.getroot().find(Tag).find('pval').find('MandatoryFieldEntries').text
            CondList = CondList.split(',')
            self.CheckPVals(DCsDict=DCsDict,
                            Log=Log,
                            FmtLog=FmtLog,
                            HeadingSpaces='    ',
                            ConditionList=CondList)
            CondList = XmlObj.getroot().find(Tag).find('HWE_pval').find('MandatoryFieldEntries').text
            CondList = CondList.split(',')
            self.CheckHWEPVals(DCsDict=DCsDict,
                               Log=Log,
                               FmtLog=FmtLog,
                               HeadingSpaces='    ',
                               ConditionList=CondList)
            CondList = XmlObj.getroot().find(Tag).find('n_total').find('MandatoryFieldEntries').text
            CondList = CondList.split(',')
            self.CheckNTotals(DCsDict=DCsDict,
                              Log=Log,
                              FmtLog=FmtLog,
                              HeadingSpaces='    ',
                              ConditionList=CondList)
            CondList = XmlObj.getroot().find(Tag).find('oevar_imp').find('MandatoryFieldEntries').text
            CondList = CondList.split(',')
            self.CheckOeVarImps(DCsDict=DCsDict,
                                Log=Log,
                                FmtLog=FmtLog,
                                HeadingSpaces='    ',
                                ConditionList=CondList)

        LogString = FmtLog.GetEndLogString()
        FmtLog.Write(LogString+'\n')
        FmtLog.Close()

        LogString  = HeadingSpaces
        LogString += '  -- Done ...'
        print LogString
        Log.Write(LogString+'\n')

        return

    def CheckUsedForImps(self,
                         DCsDict={},
                         Log=Logger,
                         FmtLog=Logger,
                         HeadingSpaces='',
                         ConditionList=[]):
        for FName, DCs in DCsDict.iteritems():
            LogString  = HeadingSpaces
            LogString += '    ++ Checking \"used_for_imp\" fields format for \"'+FName+'\" ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

            FilterArray    = (DCs.DataContainers['used_for_imp'].GetDataArray()==ConditionList[0])
            for i in range(1,len(ConditionList)):
                TmpFilterArray = (DCs.DataContainers['used_for_imp'].GetDataArray()==ConditionList[i])
                FilterArray    = (FilterArray | TmpFilterArray)

            LogString  = HeadingSpaces
            LogString += '    ** Array \"used_for_imp\" complies with the ENGAGE analysis plan v3.0!'
            boComplies = True
            if(len(scipy.compress(FilterArray,DCs.DataContainers['imputed'].GetDataArray()))!=
               len(DCs.DataContainers['imputed'].GetDataArray())):
                LogString  = HeadingSpaces
                LogString += '    ** Array \"used_for_imp\" does not comply with the ENGAGE analysis plan v3.0!'
                boComplies = False
            self.SetboFieldFormatOKDict(ColumnID='used_for_imp',
                                        Flag=boComplies)

            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

            LogString  = HeadingSpaces
            LogString += '    -- Done ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

        return

    def CheckImputeds(self,
                      DCsDict={},
                      Log=Logger,
                      FmtLog=Logger,
                      HeadingSpaces='',
                      ConditionList=[]):
        for FName, DCs in DCsDict.iteritems():
            LogString  = HeadingSpaces
            LogString += '    ++ Checking \"imputed\" fields format for \"'+FName+'\" ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

            FilterArray    = (DCs.DataContainers['imputed'].GetDataArray()==ConditionList[0])
            for i in range(1,len(ConditionList)):
                TmpFilterArray = (DCs.DataContainers['imputed'].GetDataArray()==ConditionList[i])
                FilterArray    = (FilterArray | TmpFilterArray)

            LogString  = HeadingSpaces
            LogString += '    ** Array \"imputed\" complies with the ENGAGE analysis plan v3.0!'
            boComplies = True
            if(len(scipy.compress(FilterArray,DCs.DataContainers['imputed'].GetDataArray()))!=
               len(DCs.DataContainers['imputed'].GetDataArray())):
                LogString  = HeadingSpaces
                LogString += '    ** Array \"imputed\" does not comply with the ENGAGE analysis plan v3.0!'
                boComplies = False
            self.SetboFieldFormatOKDict(ColumnID='imputed',
                                        Flag=boComplies)

            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

            LogString  = HeadingSpaces
            LogString += '    -- Done ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

        return

    def CheckStrandGenomes(self,
                           DCsDict={},
                           Log=Logger,
                           FmtLog=Logger,
                           HeadingSpaces='',
                           ConditionList=[]):
        for FName, DCs in DCsDict.iteritems():
            LogString  = HeadingSpaces
            LogString += '    ++ Checking \"strand_genome\" fields format for \"'+FName+'\" ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

            FilterArray    = (DCs.DataContainers['strand_genome'].GetDataArray()==ConditionList[0])
            for i in range(1,len(ConditionList)):
                TmpFilterArray = (DCs.DataContainers['strand_genome'].GetDataArray()==ConditionList[i])
                FilterArray    = (FilterArray | TmpFilterArray)

            LogString  = HeadingSpaces
            LogString += '    ** Array \"strand_genome\" complies with the ENGAGE analysis plan v3.0!'
            boComplies = True
            if(len(scipy.compress(FilterArray,DCs.DataContainers['strand_genome'].GetDataArray()))!=
               len(DCs.DataContainers['strand_genome'].GetDataArray())):
                LogString  = HeadingSpaces
                LogString += '    ** Array \"strand_genome\" does not comply with the ENGAGE analysis plan v3.0!'
                boComplies = False
            self.SetboFieldFormatOKDict(ColumnID='strand_genome',
                                        Flag=boComplies)

            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

            LogString  = HeadingSpaces
            LogString += '    -- Done ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

        return

    def CheckChrs(self,
                  DCsDict={},
                  Log=Logger,
                  FmtLog=Logger,
                  HeadingSpaces='',
                  ConditionList=[]):
        for FName, DCs in DCsDict.iteritems():
            LogString  = HeadingSpaces
            LogString += '    ++ Checking \"chr\" fields format for \"'+FName+'\" ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

            FilterArray    = None
            TmpFilterArray = None
            for Condition in ConditionList:
                if(re.match('RANGE',Condition)):
                    Start = int(Condition.split('[')[-1].split(']')[0].split('-')[0])
                    End   = int(Condition.split('[')[-1].split(']')[0].split('-')[-1])
                    Range = []
                    for i in range(Start,End+1):
                        Range.append(str(i))
                    if(FilterArray==None):
                        FilterArray = (DCs.DataContainers['chr'].GetDataArray()==Range[0])
                        for i  in range(1,len(Range)):
                            TmpFilterArray = (DCs.DataContainers['chr'].GetDataArray()==Range[i])
                            FilterArray    = (FilterArray | TmpFilterArray)
                    else:
                        for i  in range(0,len(Range)):
                            TmpFilterArray = (DCs.DataContainers['chr'].GetDataArray()==Range[i])
                            FilterArray    = (FilterArray | TmpFilterArray)
                else:
                    if(FilterArray==None):
                        FilterArray = DCs.DataContainers['chr'].GetDataArray()==Condition
                    else:
                        TmpFilterArray = DCs.DataContainers['chr'].GetDataArray()==Condition
                        FilterArray    = (FilterArray | TmpFilterArray)

            LogString  = HeadingSpaces
            LogString += '    ** Array \"chr\" complies with the ENGAGE analysis plan v3.0!'
            boComplies = True
            if(len(scipy.compress(FilterArray,DCs.DataContainers['chr'].GetDataArray()))!=
               len(DCs.DataContainers['chr'].GetDataArray())):
                LogString  = HeadingSpaces
                LogString += '    ** Array \"chr\" does not comply with the ENGAGE analysis plan v3.0!'
                boComplies = False
            self.SetboFieldFormatOKDict(ColumnID='chr',
                                        Flag=boComplies)

            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')
            LogString  = HeadingSpaces
            LogString += '    -- Done ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

        return

    def CheckSNPIDs(self,
                    DCsDict={},
                    Log=Logger,
                    FmtLog=Logger,
                    HeadingSpaces='',
                    ConditionList=[]):
        for FName, DCs in DCsDict.iteritems():
            LogString  = HeadingSpaces
            LogString += '    ++ Checking \"SNPID\" fields format for \"'+FName+'\" ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')
            LogString  = HeadingSpaces
            boComplies = True
            for Entry in DCs.DataContainers['SNPID'].GetDataArray():
                boConditionHolds = bool(re.search(ConditionList[0],Entry))
                for i in range(1,len(ConditionList)):
                    boConditionHolds = (boConditionHolds or bool(re.search(ConditionList[i],Entry)))
                if(not boConditionHolds):
                    LogString += '    ** SNPID does not comply with the ENGAGE analysis plan v3.0!\n'
                    boComplies = False
                    break

            self.SetboFieldFormatOKDict(ColumnID='SNPID',
                                        Flag=boComplies)
            if(boComplies):
                LogString += '    ** Array \"SNPID\" complies with the ENGAGE analysis plan v3.0!'

            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')
            LogString  = HeadingSpaces
            LogString += '    -- Done ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

        return

    def CheckColumnFormat(self,
                          DCsDict={},
                          Log=Logger,
                          FmtLog=Logger,
                          HeadingSpaces=''):
        for FName, DCs in DCsDict.iteritems():
            LogString  = HeadingSpaces
            LogString += '    ++ Checking column format for \"'+FName+'\" ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')
            for Key in DCs.DataContainers.iterkeys():
                LogString  = HeadingSpaces
                LogString += '    ** Column name \"'+Key+'\" '
                if((Key in self.GetColumnFormat()) and
                   (DCs.DataContainers[Key].GetDataName() in self.GetColumnFormat())):
                    LogString += 'complies with the ENGAGE analysis plan v3.0!'
                    self.SetboColumnFormatOK(Flag=True)
                else:
                    LogString += 'does not comply with the ENGAGE analysis plan v3.0!'
                    self.SetboColumnFormatOK(Flag=False)
                print LogString
                Log.Write(LogString+'\n')
                FmtLog.Write(LogString+'\n')
            LogString  = HeadingSpaces
            LogString += '    -- Done ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

        return

    def CheckIfFilesExist(self,
                          XmlObj=lxml.etree._ElementTree,
                          Tag=str,
                          Log=Logger,
                          HeadingSpaces=''):
        LogString  = HeadingSpaces
        LogString += '  ++ Checking if \"'+Tag+'\" exist(s) and are (is a) regular file(s) ...'
        print LogString
        Log.Write(LogString+'\n')
        for F in XmlObj.getroot().find(Tag):
            if(eval(F.find('boUse').text)):
                FName = os.path.join(F.find('Path').text,F.find('Name').text)
                if(os.path.exists(FName)):
                    LogString  = HeadingSpaces
                    LogString += '  ** Path \"'+FName+'\" exists ...'
                    print LogString
                    Log.Write(LogString+'\n')
                else:
                    LogString  = HeadingSpaces
                    LogString += '  ** Path \"'+FName+'\" does not exist ...\n'
                    LogString += '!!! EXITING !!!'
                    print LogString
                    Log.Write(LogString+'\n')
                    sys.exit(1)
                if(os.path.isfile(FName)):
                    LogString  = HeadingSpaces
                    LogString += '  ** Path \"'+FName+'\" is a regular file ...'
                    print LogString
                    Log.Write(LogString+'\n')
                else:
                    LogString  = HeadingSpaces
                    LogString += '  ** Path \"'+FName+'\" is not a regular file ...\n'
                    LogString += '!!! EXITING !!!'
                    print LogString
                    Log.Write(LogString+'\n')
                    sys.exit(1)

        LogString  = HeadingSpaces
        LogString += '  -- Done ...'
        print LogString
        Log.Write(LogString+'\n')
        return
