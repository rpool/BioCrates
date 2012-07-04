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
            self.SetExtraInfoFiles()
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
            LogString += '    ** Removed '+str(NLinesInFile-NLinesInArray)+' duplicate lines!'
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
                            Log=Logger,
                            boRemoveDuplicateLines=False):
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

            if(boRemoveDuplicateLines):
                NLinesInFile,\
                NLinesInArray  = EIFile.ParseToLineArray()
                LogString = '      ** Removed '+str(NLinesInFile-NLinesInArray)+' duplicate lines!'
                print LogString
                Log.Write(LogString+'\n')
                DCsDict[EIFile.GetName()] = EIFile.LineArray2DataContainers()
            else:
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

        for Key in DCsDict.iterkeys():
            for Column in XmlObj.getroot().find(Tag):
                if(type(Column)!=lxml.etree._Comment and
                   eval(Column.find('boSetColumnName').text)):
                    Source = Column.find('ColumnName').text
                    Dest   = Column.find('SetColumnName').text

                    LogString  = '      ++ Renaming column \"'+Source+\
                                     '\" to \"'+Dest+'\" ...'
                    print LogString
                    Log.Write(LogString+'\n')
                    FmtLog.Write(LogString+'\n')

                    DCsDict[Key].DataContainers[Dest] = DCsDict[Key].DataContainers[Source]
                    DCsDict[Key].DataContainers[Dest].RenameColumnOfDataArray(Source,
                                                                              Dest)
                    del DCsDict[Key].DataContainers[Source]

                    LogString  = '      -- Done ...'
                    print LogString
                    Log.Write(LogString+'\n')
                    FmtLog.Write(LogString+'\n')

        for Key in DCsDict.iterkeys():
            for Column in XmlObj.getroot().find(Tag):
                if(Column.find('Renames')!=None):
                    for Rename in Column.find('Renames'):
                        Source = Rename.find('Source').text
                        Dest   = Rename.find('Dest').text

                        LogString  = '      ++ Renaming \"'+Source+\
                                     '\" fields to \"'+Dest+\
                                     '\" for column \"'+Column.tag+\
                                     '\" ...'
                        print LogString
                        Log.Write(LogString+'\n')
                        FmtLog.Write(LogString+'\n')

                        DCsDict[Key].DataContainers[Column.tag].RenameFieldsInDataArray(Source,
                                                                                        Dest)

                        LogString  = '      -- Done ...'
                        print LogString
                        Log.Write(LogString+'\n')
                        FmtLog.Write(LogString+'\n')

        for Key in DCsDict.iterkeys():
            SNPIDColumn = XmlObj.getroot().find(Tag).find('SNPID').tag

            LogString  = '      ++ Checking on duplicate SNPIDs ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

            DuplicateDict = DCsDict[Key].DataContainers[SNPIDColumn].FindDuplicates()

            LogString  = '        ** The maximum number of duplicate SNPs is '+\
                         str(DCsDict[Key].DataContainers[SNPIDColumn].GetMaxNDuplicates())
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')
            if(len(DuplicateDict)>0):
                LogString  = '        ** Found the following duplicate SNPs:\n'
                for Key, Value in DuplicateDict.iteritems():
                    LogString += '           '+Key+' (Occurrence: '+str(Value)+')\n'
                print LogString[:-1]
                Log.Write(LogString[:-1]+'\n')
                FmtLog.Write(LogString[-1]+'\n')
            else:
                LogString  = '        ** No duplicates found!'
                print LogString
                Log.Write(LogString+'\n')
                FmtLog.Write(LogString+'\n')

            LogString  = '      -- Done ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

        self.CheckColumnFormat(DCsDict=DCsDict,
                               Log=Log,
                               FmtLog=FmtLog,
                               HeadingSpaces='  ')
        CondList = XmlObj.getroot().find(Tag).find('SNPID').find('MandatoryFieldEntries').text
        CondList = CondList.split(',')
        self.CheckSNPIDs(DCsDict=DCsDict,
                         Log=Log,
                         FmtLog=FmtLog,
                         HeadingSpaces='  ',
                         ConditionList=CondList)
        CondList = XmlObj.getroot().find(Tag).find('chr').find('MandatoryFieldEntries').text
        CondList = CondList.split(',')
        self.CheckChrs(DCsDict=DCsDict,
                       Log=Log,
                       FmtLog=FmtLog,
                       HeadingSpaces='  ',
                       ConditionList=CondList)
        CondList = XmlObj.getroot().find(Tag).find('strand_genome').find('MandatoryFieldEntries').text
        CondList = CondList.split(',')
        self.CheckStrandGenomes(DCsDict=DCsDict,
                                Log=Log,
                                FmtLog=FmtLog,
                                HeadingSpaces='  ',
                                ConditionList=CondList)
        CondList = XmlObj.getroot().find(Tag).find('imputed').find('MandatoryFieldEntries').text
        CondList = CondList.split(',')
        self.CheckImputeds(DCsDict=DCsDict,
                           Log=Log,
                           FmtLog=FmtLog,
                           HeadingSpaces='  ',
                           ConditionList=CondList)
        CondList = XmlObj.getroot().find(Tag).find('used_for_imp').find('MandatoryFieldEntries').text
        CondList = CondList.split(',')
        self.CheckUsedForImps(DCsDict=DCsDict,
                              Log=Log,
                              FmtLog=FmtLog,
                              HeadingSpaces='  ',
                              ConditionList=CondList)

        if(Tag=='MtbGWAColumns'):
            CondList = XmlObj.getroot().find(Tag).find('position').find('MandatoryFieldEntries').text
            CondList = CondList.split(',')
            self.CheckPositions(DCsDict=DCsDict,
                                Log=Log,
                                FmtLog=FmtLog,
                                HeadingSpaces='  ',
                                ConditionList=CondList)
            CondList = XmlObj.getroot().find(Tag).find('coded_all').find('MandatoryFieldEntries').text
            CondList = CondList.split(',')
            self.CheckCodedAlls(DCsDict=DCsDict,
                                Log=Log,
                                FmtLog=FmtLog,
                                HeadingSpaces='  ',
                                ConditionList=CondList)
            CondList = XmlObj.getroot().find(Tag).find('noncoded_all').find('MandatoryFieldEntries').text
            CondList = CondList.split(',')
            self.CheckNonCodedAlls(DCsDict=DCsDict,
                                   Log=Log,
                                   FmtLog=FmtLog,
                                   HeadingSpaces='  ',
                                   ConditionList=CondList)
            CondList   = XmlObj.getroot().find(Tag).find('beta').find('MandatoryFieldEntries').text
            CondList   = CondList.split(',')
            Precision  = int(XmlObj.getroot().find(Tag).find('beta').find('NumericFieldPrecision').text)
            APrecision = int(XmlObj.getroot().find(Tag).find('beta').find('AllowedNumericFieldPrecision').text)
            self.CheckBetas(DCsDict=DCsDict,
                            Log=Log,
                            FmtLog=FmtLog,
                            HeadingSpaces='  ',
                            ConditionList=CondList,
                            Precision=Precision,
                            AllowedPrecision=APrecision)
            CondList   = XmlObj.getroot().find(Tag).find('SE').find('MandatoryFieldEntries').text
            CondList   = CondList.split(',')
            Precision  = int(XmlObj.getroot().find(Tag).find('SE').find('NumericFieldPrecision').text)
            APrecision = int(XmlObj.getroot().find(Tag).find('SE').find('AllowedNumericFieldPrecision').text)
            self.CheckSEs(DCsDict=DCsDict,
                          Log=Log,
                          FmtLog=FmtLog,
                          HeadingSpaces='  ',
                          ConditionList=CondList,
                          Precision=Precision,
                          AllowedPrecision=APrecision)
            CondList = XmlObj.getroot().find(Tag).find('pval').find('MandatoryFieldEntries').text
            CondList = CondList.split(',')
            self.CheckPVals(DCsDict=DCsDict,
                            Log=Log,
                            FmtLog=FmtLog,
                            HeadingSpaces='  ',
                            ConditionList=CondList)
            CondList = XmlObj.getroot().find(Tag).find('HWE_pval').find('MandatoryFieldEntries').text
            CondList = CondList.split(',')
            self.CheckHWEPVals(DCsDict=DCsDict,
                               Log=Log,
                               FmtLog=FmtLog,
                               HeadingSpaces='  ',
                               ConditionList=CondList)
            CondList = XmlObj.getroot().find(Tag).find('n_total').find('MandatoryFieldEntries').text
            CondList = CondList.split(',')
            self.CheckNTotals(DCsDict=DCsDict,
                              Log=Log,
                              FmtLog=FmtLog,
                              HeadingSpaces='  ',
                              ConditionList=CondList)
            CondList = XmlObj.getroot().find(Tag).find('oevar_imp').find('MandatoryFieldEntries').text
            CondList = CondList.split(',')
            self.CheckOeVarImps(DCsDict=DCsDict,
                                Log=Log,
                                FmtLog=FmtLog,
                                HeadingSpaces='  ',
                                ConditionList=CondList)
            CondList = XmlObj.getroot().find(Tag).find('AF_coded_all').find('MandatoryFieldEntries').text
            CondList = CondList.split(',')
            self.CheckAfCodedAlls(DCsDict=DCsDict,
                                  Log=Log,
                                  FmtLog=FmtLog,
                                  HeadingSpaces='  ',
                                  ConditionList=CondList)

        LogString = FmtLog.GetEndLogString()
        FmtLog.Write(LogString+'\n')
        FmtLog.Close()

        LogString  = HeadingSpaces
        LogString += '  -- Done ...'
        print LogString
        Log.Write(LogString+'\n')

        return

    def CheckAfCodedAlls(self,
                         DCsDict={},
                         Log=Logger,
                         FmtLog=Logger,
                         HeadingSpaces='',
                         ConditionList=[]):

        for FName, DCs in DCsDict.iteritems():
            LogString  = HeadingSpaces
            LogString += '    ++ Checking \"AF_coded_all\" fields format for \"'+FName+'\" ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

            FilterArray = None
            for i in range(len(ConditionList)):
                if((ConditionList[i][0]=='[') and
                   (ConditionList[i][-1]==']')):
                    RegExp = re.compile(ConditionList[i])
                    VMatch = scipy.vectorize(lambda x:bool(RegExp.match(x)))
                    if(i==0):
                        FilterArray = VMatch(DCs.DataContainers['AF_coded_all'].GetDataArray())
                    else:
                        TmpFilterArray = VMatch(DCs.DataContainers['AF_coded_all'].GetDataArray())
                        FilterArray    = (FilterArray | TmpFilterArray)
                else:
                    if(i==0):
                        FilterArray = (DCs.DataContainers['AF_coded_all'].GetDataArray()==ConditionList[i])
                    else:
                        TmpFilterArray = (DCs.DataContainers['AF_coded_all'].GetDataArray()==ConditionList[i])
                        FilterArray    = (FilterArray | TmpFilterArray)

            LogString  = HeadingSpaces
            LogString += '      ** Array \"AF_coded_all\" complies with the ENGAGE analysis plan v3.0!'
            boComplies = True
            if(len(scipy.compress(FilterArray,DCs.DataContainers['AF_coded_all'].GetDataArray()))!=
               len(DCs.DataContainers['AF_coded_all'].GetDataArray())):
                LogString  = HeadingSpaces
                LogString += '      ** Array \"AF_coded_all\" does not comply with the ENGAGE analysis plan v3.0!'
                boComplies = False
            self.SetboFieldFormatOKDict(ColumnID='AF_coded_all',
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

    def CheckOeVarImps(self,
                       DCsDict={},
                       Log=Logger,
                       FmtLog=Logger,
                       HeadingSpaces='',
                       ConditionList=[]):

        for FName, DCs in DCsDict.iteritems():
            LogString  = HeadingSpaces
            LogString += '    ++ Checking \"oevar_imp\" fields format for \"'+FName+'\" ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

            FilterArray = None
            for i in range(len(ConditionList)):
                if((ConditionList[i][0]=='[') and
                   (ConditionList[i][-1]==']')):
                    RegExp = re.compile(ConditionList[i])
                    VSearch = scipy.vectorize(lambda x:bool(RegExp.search(x)))
                    if(i==0):
                        FilterArray = VSearch(DCs.DataContainers['oevar_imp'].GetDataArray())
                    else:
                        TmpFilterArray = VSearch(DCs.DataContainers['oevar_imp'].GetDataArray())
                        FilterArray    = (FilterArray | TmpFilterArray)
                else:
                    if(i==0):
                        FilterArray = (DCs.DataContainers['oevar_imp'].GetDataArray()==ConditionList[i])
                    else:
                        TmpFilterArray = (DCs.DataContainers['oevar_imp'].GetDataArray()==ConditionList[i])
                        FilterArray    = (FilterArray | TmpFilterArray)

            LogString  = HeadingSpaces
            LogString += '      ** Array \"oevar_imp\" complies with the ENGAGE analysis plan v3.0!'
            boComplies = True
            if(len(scipy.compress(FilterArray,DCs.DataContainers['oevar_imp'].GetDataArray()))!=
               len(DCs.DataContainers['oevar_imp'].GetDataArray())):
                LogString  = HeadingSpaces
                LogString += '      ** Array \"oevar_imp\" does not comply with the ENGAGE analysis plan v3.0!'
                boComplies = False
            self.SetboFieldFormatOKDict(ColumnID='oevar_imp',
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

    def CheckNTotals(self,
                     DCsDict={},
                     Log=Logger,
                     FmtLog=Logger,
                     HeadingSpaces='',
                     ConditionList=[]):

        for FName, DCs in DCsDict.iteritems():
            LogString  = HeadingSpaces
            LogString += '    ++ Checking \"n_total\" fields format for \"'+FName+'\" ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

            FilterArray = None
            for i in range(len(ConditionList)):
                if((ConditionList[i][0]=='[') and
                   (ConditionList[i][-1]==']')):
                    RegExp = re.compile(ConditionList[i])
                    VMatch = scipy.vectorize(lambda x:bool(RegExp.match(x)))
                    if(i==0):
                        FilterArray = VMatch(DCs.DataContainers['n_total'].GetDataArray())
                    else:
                        TmpFilterArray = VMatch(DCs.DataContainers['n_total'].GetDataArray())
                        FilterArray    = (FilterArray | TmpFilterArray)
                else:
                    if(i==0):
                        FilterArray = (DCs.DataContainers['n_total'].GetDataArray()==ConditionList[i])
                    else:
                        TmpFilterArray = (DCs.DataContainers['n_total'].GetDataArray()==ConditionList[i])
                        FilterArray    = (FilterArray | TmpFilterArray)

            LogString  = HeadingSpaces
            LogString += '      ** Array \"n_total\" complies with the ENGAGE analysis plan v3.0!'
            boComplies = True
            if(len(scipy.compress(FilterArray,DCs.DataContainers['n_total'].GetDataArray()))!=
               len(DCs.DataContainers['n_total'].GetDataArray())):
                LogString  = HeadingSpaces
                LogString += '      ** Array \"n_total\" does not comply with the ENGAGE analysis plan v3.0!'
                boComplies = False
            self.SetboFieldFormatOKDict(ColumnID='n_total',
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

    def CheckPVals(self,
                   DCsDict={},
                   Log=Logger,
                   FmtLog=Logger,
                   HeadingSpaces='',
                   ConditionList=[]):

        for FName, DCs in DCsDict.iteritems():
            LogString  = HeadingSpaces
            LogString += '    ++ Checking \"pval\" fields format for \"'+FName+'\" ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

            FilterArray = None
            for i in range(len(ConditionList)):
                RegExp  = re.compile(ConditionList[i])
                VSearch = scipy.vectorize(lambda x:bool(RegExp.search(x)))
                if(i==0):
                    FilterArray = VSearch(DCs.DataContainers['pval'].GetDataArray())
                else:
                    TmpFilterArray = VSearch(DCs.DataContainers['pval'].GetDataArray())
                    FilterArray    = (FilterArray | TmpFilterArray)

            LogString  = HeadingSpaces
            LogString += '      ** Array \"pval\" complies with the ENGAGE analysis plan v3.0!'
            boComplies = True
            if(len(scipy.compress(FilterArray,DCs.DataContainers['pval'].GetDataArray()))!=
               len(DCs.DataContainers['pval'].GetDataArray())):
                LogString  = HeadingSpaces
                LogString += '      ** Array \"pval\" does not comply with the ENGAGE analysis plan v3.0!'
                boComplies = False
            self.SetboFieldFormatOKDict(ColumnID='pval',
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

    def CheckHWEPVals(self,
                      DCsDict={},
                      Log=Logger,
                      FmtLog=Logger,
                      HeadingSpaces='',
                      ConditionList=[]):

        for FName, DCs in DCsDict.iteritems():
            LogString  = HeadingSpaces
            LogString += '    ++ Checking \"HWE_pval\" fields format for \"'+FName+'\" ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

            FilterArray = None
            for i in range(len(ConditionList)):
                RegExp  = re.compile(ConditionList[i])
                VSearch = scipy.vectorize(lambda x:bool(RegExp.search(x)))
                if(i==0):
                    FilterArray = VSearch(DCs.DataContainers['HWE_pval'].GetDataArray())
                else:
                    TmpFilterArray = VSearch(DCs.DataContainers['HWE_pval'].GetDataArray())
                    FilterArray    = (FilterArray | TmpFilterArray)

            LogString  = HeadingSpaces
            LogString += '      ** Array \"HWE_pval\" complies with the ENGAGE analysis plan v3.0!'
            boComplies = True
            if(len(scipy.compress(FilterArray,DCs.DataContainers['HWE_pval'].GetDataArray()))!=
               len(DCs.DataContainers['HWE_pval'].GetDataArray())):
                LogString  = HeadingSpaces
                LogString += '      ** Array \"HWE_pval\" does not comply with the ENGAGE analysis plan v3.0!'
                boComplies = False
            self.SetboFieldFormatOKDict(ColumnID='HWE_pval',
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

    def CheckSEs(self,
                 DCsDict={},
                 Log=Logger,
                 FmtLog=Logger,
                 HeadingSpaces='',
                 ConditionList=[],
                 Precision=int,
                 AllowedPrecision=int):

        for FName, DCs in DCsDict.iteritems():
            LogString  = HeadingSpaces
            LogString += '    ++ Checking \"SE\" fields format for \"'+FName+'\" ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

            FilterArray = None
            for i in range(len(ConditionList)):
                RegExp  = re.compile(ConditionList[i])
                VSearch = scipy.vectorize(lambda x:bool(RegExp.search(x)))
                if(i==0):
                    FilterArray = VSearch(DCs.DataContainers['SE'].GetDataArray())
                else:
                    TmpFilterArray = VSearch(DCs.DataContainers['SE'].GetDataArray())
                    FilterArray    = (FilterArray | TmpFilterArray)
            TmpFilterArray = []
            boPrecisionOK  = True
            MinPrecision   = 500
            for Entry in DCs.DataContainers['SE'].GetDataArray():
                ESplit = Entry.split('.')[-1]
                ESplit = ESplit.split('e')[0]
                ESplit = ESplit.split('E')[0]
                if(len(ESplit)>=Precision):
                    TmpFilterArray.append(True)
                else:
                    TmpFilterArray.append(False)
                    boPrecisionOK = False
                MinPrecision = min(MinPrecision,len(ESplit))
            if(MinPrecision>AllowedPrecision):
                TmpFilterArray = scipy.array(TmpFilterArray)
                FilterArray    = (FilterArray & TmpFilterArray)
            LogString = HeadingSpaces
            if(boPrecisionOK):
                LogString += '      ** Array \"SE\" complies with the precision prerequisite of the ENGAGE analysis plan v3.0!'
            else:
                LogString += '      ** Array \"SE\" does not comply with the precision prerequisite of the ENGAGE analysis plan v3.0!'
            LogString += '\n'
            LogString += '        ** Minimal precision is: '+str(MinPrecision)+' ...\n'
            LogString += '        ** Precision override value is: '+str(AllowedPrecision)+' ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

            LogString  = HeadingSpaces
            LogString += '      ** Array \"SE\" complies with the ENGAGE analysis plan v3.0!'
            boComplies = True
            if(len(scipy.compress(FilterArray,DCs.DataContainers['SE'].GetDataArray()))!=
               len(DCs.DataContainers['SE'].GetDataArray())):
                LogString  = HeadingSpaces
                LogString += '      ** Array \"SE\" does not comply with the ENGAGE analysis plan v3.0!'
                boComplies = False
            self.SetboFieldFormatOKDict(ColumnID='SE',
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

    def CheckBetas(self,
                   DCsDict={},
                   Log=Logger,
                   FmtLog=Logger,
                   HeadingSpaces='',
                   ConditionList=[],
                   Precision=int,
                   AllowedPrecision=int):

        for FName, DCs in DCsDict.iteritems():
            LogString  = HeadingSpaces
            LogString += '    ++ Checking \"beta\" fields format for \"'+FName+'\" ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

            FilterArray = None
            for i in range(len(ConditionList)):
                RegExp  = re.compile(ConditionList[i])
                VSearch = scipy.vectorize(lambda x:bool(RegExp.search(x)))
                if(i==0):
                    FilterArray = VSearch(DCs.DataContainers['beta'].GetDataArray())
                else:
                    TmpFilterArray = VSearch(DCs.DataContainers['beta'].GetDataArray())
                    FilterArray    = (FilterArray | TmpFilterArray)

            TmpFilterArray = []
            boPrecisionOK  = True
            MinPrecision   = 500
            for Entry in DCs.DataContainers['beta'].GetDataArray():
                ESplit = Entry.split('.')[-1]
                ESplit = ESplit.split('e')[0]
                ESplit = ESplit.split('E')[0]
                if(len(ESplit)>=5):
                    TmpFilterArray.append(True)
                else:
                    TmpFilterArray.append(False)
                    boPrecisionOK = False
                MinPrecision = min(MinPrecision,len(ESplit))
            if(MinPrecision>AllowedPrecision):
                TmpFilterArray = scipy.array(TmpFilterArray)
                FilterArray    = (FilterArray & TmpFilterArray)
            LogString = HeadingSpaces
            if(boPrecisionOK):
                LogString += '      ** Array \"beta\" complies with the precision prerequisite of the ENGAGE analysis plan v3.0!'
            else:
                LogString += '      ** Array \"beta\" does not comply with the precision prerequisite of the ENGAGE analysis plan v3.0!'
            LogString += '\n'
            LogString += '        ** Minimal precision is: '+str(MinPrecision)+' ...\n'
            LogString += '        ** Precision override value is: '+str(AllowedPrecision)+' ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

            LogString  = HeadingSpaces
            LogString += '      ** Array \"beta\" complies with the ENGAGE analysis plan v3.0!'
            boComplies = True
            if(len(scipy.compress(FilterArray,DCs.DataContainers['beta'].GetDataArray()))!=
               len(DCs.DataContainers['beta'].GetDataArray())):
                LogString  = HeadingSpaces
                LogString += '      ** Array \"beta\" does not comply with the ENGAGE analysis plan v3.0!'
                boComplies = False
            self.SetboFieldFormatOKDict(ColumnID='beta',
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

    def CheckNonCodedAlls(self,
                          DCsDict={},
                          Log=Logger,
                          FmtLog=Logger,
                          HeadingSpaces='',
                          ConditionList=[]):

        for FName, DCs in DCsDict.iteritems():
            LogString  = HeadingSpaces
            LogString += '    ++ Checking \"noncoded_all\" fields format for \"'+FName+'\" ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

            FilterArray    = (DCs.DataContainers['noncoded_all'].GetDataArray()==ConditionList[0])
            for i in range(1,len(ConditionList)):
                TmpFilterArray = (DCs.DataContainers['noncoded_all'].GetDataArray()==ConditionList[i])
                FilterArray    = (FilterArray | TmpFilterArray)

            LogString  = HeadingSpaces
            LogString += '      ** Array \"noncoded_all\" complies with the ENGAGE analysis plan v3.0!'
            boComplies = True
            if(len(scipy.compress(FilterArray,DCs.DataContainers['noncoded_all'].GetDataArray()))!=
               len(DCs.DataContainers['noncoded_all'].GetDataArray())):
                LogString  = HeadingSpaces
                LogString += '      ** Array \"noncoded_all\" does not comply with the ENGAGE analysis plan v3.0!'
                boComplies = False
            self.SetboFieldFormatOKDict(ColumnID='noncoded_all',
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

    def CheckCodedAlls(self,
                       DCsDict={},
                       Log=Logger,
                       FmtLog=Logger,
                       HeadingSpaces='',
                       ConditionList=[]):

        for FName, DCs in DCsDict.iteritems():
            LogString  = HeadingSpaces
            LogString += '    ++ Checking \"coded_all\" fields format for \"'+FName+'\" ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

            FilterArray    = (DCs.DataContainers['coded_all'].GetDataArray()==ConditionList[0])
            for i in range(1,len(ConditionList)):
                TmpFilterArray = (DCs.DataContainers['coded_all'].GetDataArray()==ConditionList[i])
                FilterArray    = (FilterArray | TmpFilterArray)

            LogString  = HeadingSpaces
            LogString += '      ** Array \"coded_all\" complies with the ENGAGE analysis plan v3.0!'
            boComplies = True
            if(len(scipy.compress(FilterArray,DCs.DataContainers['coded_all'].GetDataArray()))!=
               len(DCs.DataContainers['coded_all'].GetDataArray())):
                LogString  = HeadingSpaces
                LogString += '      ** Array \"coded_all\" does not comply with the ENGAGE analysis plan v3.0!'
                boComplies = False
            self.SetboFieldFormatOKDict(ColumnID='coded_all',
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

    def CheckPositions(self,
                       DCsDict={},
                       Log=Logger,
                       FmtLog=Logger,
                       HeadingSpaces='',
                       ConditionList=[]):

        for FName, DCs in DCsDict.iteritems():
            LogString  = HeadingSpaces
            LogString += '    ++ Checking \"position\" fields format for \"'+FName+'\" ...'
            print LogString
            Log.Write(LogString+'\n')
            FmtLog.Write(LogString+'\n')

            FilterArray = None
            for i in range(len(ConditionList)):
                if((ConditionList[i][0]=='[') and
                   (ConditionList[i][-1]==']')):
                    RegExp = re.compile(ConditionList[i])
                    VMatch = scipy.vectorize(lambda x:bool(RegExp.match(x)))
                    if(i==0):
                        FilterArray = VMatch(DCs.DataContainers['position'].GetDataArray())
                    else:
                        TmpFilterArray = VMatch(DCs.DataContainers['position'].GetDataArray())
                        FilterArray    = (FilterArray | TmpFilterArray)
                else:
                    if(i==0):
                        FilterArray = (DCs.DataContainers['position'].GetDataArray()==ConditionList[i])
                    else:
                        TmpFilterArray = (DCs.DataContainers['position'].GetDataArray()==ConditionList[i])
                        FilterArray    = (FilterArray | TmpFilterArray)

            LogString  = HeadingSpaces
            LogString += '      ** Array \"position\" complies with the ENGAGE analysis plan v3.0!'
            boComplies = True
            if(len(scipy.compress(FilterArray,DCs.DataContainers['position'].GetDataArray()))!=
               len(DCs.DataContainers['position'].GetDataArray())):
                LogString  = HeadingSpaces
                LogString += '      ** Array \"position\" does not comply with the ENGAGE analysis plan v3.0!'
                boComplies = False
            self.SetboFieldFormatOKDict(ColumnID='position',
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
            LogString += '      ** Array \"used_for_imp\" complies with the ENGAGE analysis plan v3.0!'
            boComplies = True
            if(len(scipy.compress(FilterArray,DCs.DataContainers['used_for_imp'].GetDataArray()))!=
               len(DCs.DataContainers['used_for_imp'].GetDataArray())):
                LogString  = HeadingSpaces
                LogString += '      ** Array \"used_for_imp\" does not comply with the ENGAGE analysis plan v3.0!'
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
            LogString += '      ** Array \"imputed\" complies with the ENGAGE analysis plan v3.0!'
            boComplies = True
            if(len(scipy.compress(FilterArray,DCs.DataContainers['imputed'].GetDataArray()))!=
               len(DCs.DataContainers['imputed'].GetDataArray())):
                LogString  = HeadingSpaces
                LogString += '      ** Array \"imputed\" does not comply with the ENGAGE analysis plan v3.0!'
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
            LogString += '      ** Array \"strand_genome\" complies with the ENGAGE analysis plan v3.0!'
            boComplies = True
            if(len(scipy.compress(FilterArray,DCs.DataContainers['strand_genome'].GetDataArray()))!=
               len(DCs.DataContainers['strand_genome'].GetDataArray())):
                LogString  = HeadingSpaces
                LogString += '      ** Array \"strand_genome\" does not comply with the ENGAGE analysis plan v3.0!'
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
            LogString += '      ** Array \"chr\" complies with the ENGAGE analysis plan v3.0!'
            boComplies = True
            if(len(scipy.compress(FilterArray,DCs.DataContainers['chr'].GetDataArray()))!=
               len(DCs.DataContainers['chr'].GetDataArray())):
                LogString  = HeadingSpaces
                LogString += '      ** Array \"chr\" does not comply with the ENGAGE analysis plan v3.0!'
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
                    LogString += '      ** SNPID does not comply with the ENGAGE analysis plan v3.0!'
                    print Entry
                    boComplies = False
                    break

            self.SetboFieldFormatOKDict(ColumnID='SNPID',
                                        Flag=boComplies)
            if(boComplies):
                LogString += '      ** Array \"SNPID\" complies with the ENGAGE analysis plan v3.0!'

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
                LogString += '      ** Column name \"'+Key+'\" '
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
                    LogString += '    ** Path \"'+FName+'\" exists ...'
                    print LogString
                    Log.Write(LogString+'\n')
                else:
                    LogString  = HeadingSpaces
                    LogString += '    ** Path \"'+FName+'\" does not exist ...\n'
                    LogString += '!!! EXITING !!!'
                    print LogString
                    Log.Write(LogString+'\n')
                    sys.exit(1)
                if(os.path.isfile(FName)):
                    LogString  = HeadingSpaces
                    LogString += '    ** Path \"'+FName+'\" is a regular file ...'
                    print LogString
                    Log.Write(LogString+'\n')
                else:
                    LogString  = HeadingSpaces
                    LogString += '    ** Path \"'+FName+'\" is not a regular file ...\n'
                    LogString += '!!! EXITING !!!'
                    print LogString
                    Log.Write(LogString+'\n')
                    sys.exit(1)

        LogString  = HeadingSpaces
        LogString += '  -- Done ...'
        print LogString
        Log.Write(LogString+'\n')
        return
