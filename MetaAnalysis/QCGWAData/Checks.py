import lxml.etree
import scipy
import re

class Checks:
    def __init__(self):
        # 0: False
        # 1: True
        # 3: Unknown
        self.MtbnameOK         = 0
        self.FormattingOK      = 0
        self.ChromosomesOK     = 0
        self.NtotalOK          = 3
        self.ImpQualOK         = 0
        self.NtotImpOK         = 0
        self.LambdaOK          = 0
        self.ScatterFreqsOK    = 3
        self.MaxNDuplicateSNPs = 0

        return

    def SetMaxNDuplicateSNPs(self,
                             Value=int):
        self.MaxNDuplicateSNPs = Value
        return

    def GetMaxNDuplicateSNPs(self):
        return self.MaxNDuplicateSNPs

    def SetMtbnameOK(self,
                     Value=int):
        self.MtbnameOK = Value
        return

    def GetMtbnameOK(self):
        return self.MtbnameOK

    def SetFormattingOK(self,
                        Value=int):
        self.FormattingOK = Value
        return

    def GetFormattingOK(self):
        return self.FormattingOK

    def SetChromosomesOK(self,
                         Value=int):
        self.ChromosomesOK = Value
        return

    def CheckChromosomesOK(self,
                           XmlObj=lxml.etree._ElementTree,
                           DataArray=scipy.array,
                           ColumnTag=str):
        CheckTag    = XmlObj.getroot().find('MtbGWAColumns').find('chr').find('Check').text
        CompareList = XmlObj.getroot().find('QCChecks').find(CheckTag).find('Compare').text
        CompareList = re.sub('\[','',CompareList)
        CompareList = re.sub('\]','',CompareList)
        CompareList = CompareList.split(',')

        UniqueDataList = scipy.unique(DataArray).tolist()

        boOk = True
        for Entry in CompareList:
            if(not (Entry in UniqueDataList)):
                boOk = False
                break

        if(boOk):
            self.SetChromosomesOK(1)

        return boOk

    def GetChromosomesOK(self):
        return self.ChromosomesOK

    def SetNtotalOK(self,
                    Value=int):
        self.NtotalOK = Value
        return

    def GetNtotalOK(self):
        return self.NtotalOK

    def SetImpQualOK(self,
                     Value=int):
        self.ImpQualOK = Value
        return

    def GetImpQualOK(self):
        return self.ImpQualOK

    def SetNtotImpOK(self,
                     Value=int):
        self.NtotImpOK = Value
        return

    def GetNtotImpOK(self):
        return self.NtotImpOK

    def SetLambdaOK(self,
                    Value=int):
        self.LambdaOK = Value
        return

    def GetLambdaOK(self):
        return self.LambdaOK

    def SetScatterFreqsOK(self,
                          Value=int):
        self.ScatterFreqsOK = Value
        return

    def GetScatterFreqsOK(self):
        return self.ScatterFreqsOK
