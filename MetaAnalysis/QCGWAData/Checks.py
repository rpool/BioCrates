import lxml.etree
import scipy
import scipy.stats
import re

import Plotting

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
        self.TTestOK           = 0
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
        CheckTag    = XmlObj.getroot().find('MtbGWAColumns').find(ColumnTag).find('Check').text
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

    def SetNTotalOK(self,
                    Value=int):
        self.NtotalOK = Value
        return

    def CheckNTotalOK(self,
                      XmlObj=lxml.etree._ElementTree,
                      DataArray=scipy.array,
                      CheckDataArray=scipy.array,
                      ColumnTag=str,
                      boPlot=False,
                      MtbName=str):
        CheckTag     = XmlObj.getroot().find('MtbGWAColumns').find(ColumnTag).find('Check').text
        CompareValue = XmlObj.getroot().find('QCChecks').find(CheckTag).find('Compare').text
        CompareType  = XmlObj.getroot().find('QCChecks').find(CheckTag).find('CompareType').text

        boOk = True

        FilterArray    = (DataArray!='NA')
        FilterArray   *= (CheckDataArray!='NA')
        DataArray      = scipy.compress(FilterArray,DataArray).astype(float)
        CheckDataArray = scipy.compress(FilterArray,CheckDataArray).astype(float)

        CorrCoeff = scipy.corrcoef(CheckDataArray,
                                   DataArray)
        if(CorrCoeff[0][1]<float(CompareValue)):
            boOk = False

        pCheckDataArray = -scipy.log10(CheckDataArray)
        pDataArray      = -scipy.log10(DataArray)
        Log10CorrCoeff  = scipy.corrcoef(pCheckDataArray,
                                         pDataArray)
        if(Log10CorrCoeff[0][1]<float(CompareValue)):
            boOk = False

        if(boPlot):
            Minimum  = 0.0
            Minimum  = min(Minimum,DataArray.min())
            Minimum  = min(Minimum,CheckDataArray.min())
            Maximum  = 1.0
            Maximum  = max(Maximum,DataArray.max())
            Maximum  = max(Maximum,CheckDataArray.max())
            Maximum += 0.05
            Plotting.Scatter(X=DataArray,
                             Y=CheckDataArray,
                             Min=Minimum,
                             Max=Maximum,
                             Legend=r'$r='+str(round(CorrCoeff[0][1],5))+'$',
                             Name='WaldPScatter_'+MtbName+'.png',
                             XLabel=r'$p-{\rm value} {~\rm (reported)} {~\rm [-]}$',
                             YLabel=r'$p-{\rm value} {~\rm (expected)} {~\rm [-]}$',
                             boDrawY_EQ_X=True)
            Minimum  = 0.0
            Minimum  = min(Minimum,pDataArray.min())
            Minimum  = min(Minimum,pCheckDataArray.min())
            Maximum  = 1.0
            Maximum  = max(Maximum,pDataArray.max())
            Maximum  = max(Maximum,pCheckDataArray.max())
            Maximum += 2.5
            Plotting.Scatter(X=pDataArray,
                             Y=pCheckDataArray,
                             Min=Minimum,
                             Max=Maximum,
                             Legend=r'$r='+str(round(Log10CorrCoeff[0][1],5))+'$',
                             Name='WaldpPScatter_'+MtbName+'.png',
                             XLabel=r'$-\log_{10}{(p-{\rm value})} {~\rm (reported)} {~\rm [-]}$',
                             YLabel=r'$-\log_{10}{(p-{\rm value})} {~\rm (expected)} {~\rm [-]}$',
                             boDrawY_EQ_X=True)

        if(boOk):
            self.SetNTotalOK(1)

        return boOk

    def GetNTotalOK(self):
        return self.NtotalOK

    def SetTTestOK(self,
                    Value=int):
        self.TTestOK = Value
        return

    def CheckTTestOK(self,
                     XmlObj=lxml.etree._ElementTree,
                     DataArray=scipy.array,
                     CheckDataArray=scipy.array,
                     ColumnTag=str,
                     boPlot=False,
                     MtbName=str):
        CheckTag     = XmlObj.getroot().find('MtbGWAColumns').find(ColumnTag).find('Check').text
        CompareValue = XmlObj.getroot().find('QCChecks').find(CheckTag).find('Compare').text
        CompareType  = XmlObj.getroot().find('QCChecks').find(CheckTag).find('CompareType').text

        boOk = True

        FilterArray    = (DataArray!='NA')
        FilterArray   *= (CheckDataArray!='NA')
        DataArray      = scipy.compress(FilterArray,DataArray).astype(float)
        CheckDataArray = scipy.compress(FilterArray,CheckDataArray).astype(float)

        CorrCoeff = scipy.corrcoef(CheckDataArray,
                                   DataArray)
        if(CorrCoeff[0][1]<float(CompareValue)):
            boOk = False

        pCheckDataArray = -scipy.log10(CheckDataArray)
        pDataArray      = -scipy.log10(DataArray)
        Log10CorrCoeff  = scipy.corrcoef(pCheckDataArray,
                                         pDataArray)
        if(Log10CorrCoeff[0][1]<float(CompareValue)):
            boOk = False

        if(boPlot):
            Minimum  = 0.0
            Minimum  = min(Minimum,DataArray.min())
            Minimum  = min(Minimum,CheckDataArray.min())
            Maximum  = 1.0
            Maximum  = max(Maximum,DataArray.max())
            Maximum  = max(Maximum,CheckDataArray.max())
            Maximum += 0.05
            Plotting.Scatter(X=DataArray,
                             Y=CheckDataArray,
                             Min=Minimum,
                             Max=Maximum,
                             Legend=r'$r='+str(round(CorrCoeff[0][1],5))+'$',
                             Name='TTestPScatter_'+MtbName+'.png',
                             XLabel=r'$p-{\rm value} {~\rm (reported)} {~\rm [-]}$',
                             YLabel=r'$p-{\rm value} {~\rm (expected)} {~\rm [-]}$',
                             boDrawY_EQ_X=True)
            Minimum  = 0.0
            Minimum  = min(Minimum,pDataArray.min())
            Minimum  = min(Minimum,pCheckDataArray.min())
            Maximum  = 1.0
            Maximum  = max(Maximum,pDataArray.max())
            Maximum  = max(Maximum,pCheckDataArray.max())
            Maximum += 2.5
            Plotting.Scatter(X=pDataArray,
                             Y=pCheckDataArray,
                             Min=Minimum,
                             Max=Maximum,
                             Legend=r'$r='+str(round(Log10CorrCoeff[0][1],5))+'$',
                             Name='TTestpPScatter_'+MtbName+'.png',
                             XLabel=r'$-\log_{10}{(p-{\rm value})} {~\rm (reported)} {~\rm [-]}$',
                             YLabel=r'$-\log_{10}{(p-{\rm value})} {~\rm (expected)} {~\rm [-]}$',
                             boDrawY_EQ_X=True)

        if(boOk):
            self.SetTTestOK(1)

        return boOk

    def GetTTestOK(self):
        return self.TTestOK

    def SetImpQualOK(self,
                     Value=int):
        self.ImpQualOK = Value
        return

    def CheckImpQualOK(self,
                       XmlObj=lxml.etree._ElementTree,
                       ImpQArray=scipy.array,
                       ImputedArray=scipy.array,
                       ImpQColumnTag=str,
                       ImputedColumnTag=str):
        ImpQCheckTag     = XmlObj.getroot().find('MtbGWAColumns').find(ImpQColumnTag).find('Check').text
        ImpQCompareValue = XmlObj.getroot().find('QCChecks').find(ImpQCheckTag).find('Compare').text
        ImpQCompareType  = XmlObj.getroot().find('QCChecks').find(ImpQCheckTag).find('CompareType').text

        ImpQFilterArray  = (ImpQArray==ImpQCompareValue)
        ImpQFilteredSize = len(scipy.compress(ImpQFilterArray,ImpQArray))

        ImputedCheckTag     = XmlObj.getroot().find('MtbGWAColumns').find(ImputedColumnTag).find('CheckImpQualOK').text
        ImputedCompareValue = XmlObj.getroot().find('QCChecks').find(ImputedCheckTag).find('Compare').text
        ImputedCompareType  = XmlObj.getroot().find('QCChecks').find(ImputedCheckTag).find('CompareType').text

        ImputedFilterArray  = (ImputedArray==ImputedCompareValue)
        ImputedFilteredSize = len(scipy.compress(ImputedFilterArray,ImputedArray))

        boOk = False
        if(ImpQFilteredSize==ImputedFilteredSize):
            if(len(scipy.compress(ImpQFilterArray*ImputedFilterArray,ImpQArray))==ImpQFilteredSize):
                boOk = True
                self.SetImpQualOK(1)

        return boOk

    def GetImpQualOK(self):

        return self.ImpQualOK

    def SetNTotImpOK(self,
                     Value=int):
        self.NtotImpOK = Value
        return

    def CheckNTotImpOK(self,
                       XmlObj=lxml.etree._ElementTree,
                       NTotalArray=scipy.array,
                       ImputedArray=scipy.array,
                       NTotalColumnTag=str,
                       ImputedColumnTag=str):
        ImputedCheckTag     = XmlObj.getroot().find('MtbGWAColumns').find(ImputedColumnTag).find('CheckNTotImpOK').text
        ImputedCompareValue = XmlObj.getroot().find('QCChecks').find(ImputedCheckTag).find('Compare').text
        ImputedCompareType  = XmlObj.getroot().find('QCChecks').find(ImputedCheckTag).find('CompareType').text

        ImputedFilterArray  = (ImputedArray==ImputedCompareValue)
        NTotalFilteredArray = scipy.compress(ImputedFilterArray,NTotalArray)

        boOk = False
        if(len(scipy.unique(NTotalFilteredArray))==1):
            boOk = True
            self.SetImpQualOK(1)

        return boOk

    def GetNTotImpOK(self):
        return self.NtotImpOK

    def SetLambdaOK(self,
                    Value=int):
        self.LambdaOK = Value
        return

    def CheckLambdaOK(self,
                      XmlObj=lxml.etree._ElementTree,
                      DataArray=scipy.array,
                      EMACArray=scipy.array,
                      MaxNTotal=int,
                      ColumnTag=str,
                      boPlot=False,
                      MtbName=str):
        CheckTag     = XmlObj.getroot().find('MtbGWAColumns').find(ColumnTag).find('Check').text
        CompareValue = XmlObj.getroot().find('QCChecks').find(CheckTag).find('Compare').text
        CompareType  = XmlObj.getroot().find('QCChecks').find(CheckTag).find('CompareType').text

        EMACLevels = XmlObj.getroot().find('EMAC').find('Strata').text.split(',')
        for i in range(len(EMACLevels)):
            Level = EMACLevels[i]
            if(re.search('[0-9]',Level)):
                EMACLevels[i] = float(Level)
            elif(Level=='MAX'):
                EMACLevels[i] = float(MaxNTotal)
        Colors = XmlObj.getroot().find('EMAC').find('Colors').text.split(',')

        boOk = True

        LambdaEst   = None
        SELambdaEst = None
        if(boPlot):
            LPValObsArray = -scipy.log10(DataArray)
#           The 1's are for df=1
            PValExpArray  = scipy.stats.chi2.sf(scipy.stats.chi2.rvs(1,\
                                                                      size=len(DataArray)),\
                                                1,)
            LPValExpArray = -scipy.log10(PValExpArray)
            LambdaEst,\
            SELambdaEst   = Plotting.PlotQQFilteredOnScore(MtbName=MtbName,
                                                           LPValExpArray=LPValExpArray,
                                                           LPValObsArray=LPValObsArray,
                                                           PValObsArray=DataArray,
                                                           EMACArray=EMACArray,
                                                           EMACLevels=EMACLevels,
                                                           Colors=Colors)
        else:
            LambdaEst,\
            SELambdaEst = Plotting.LambdaEstimate(Array=DataArray,
                                                  Filter=True)
        if(LambdaEst>=float(CompareValue)):
            boOk = False
        if(boOk):
            self.SetLambdaOK(1)

        return boOk,\
               LambdaEst,\
               SELambdaEst

    def GetLambdaOK(self):
        return self.LambdaOK

    def SetScatterFreqsOK(self,
                          Value=int):
        self.ScatterFreqsOK = Value
        return

    def GetScatterFreqsOK(self):
        return self.ScatterFreqsOK
