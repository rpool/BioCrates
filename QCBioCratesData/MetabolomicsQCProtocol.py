import scipy

class MetabolomicsQCProtocol:
    """
    Taken from
    Flagship Metabolomics ENGAGE - Protocol for GWAS and Meta-Analysis of Biocrates Metabolites
    Version 2.0
    """
    def __init__(self):
#       1a)
        self.RefererenceMeanCVOfAllPlatesThreshold = None
        self.MetaboliteMissingValueRateThreshold   = None
#       1b)
        self.MetaboliteOutlyingSampleThreshold     = None
        self.NOutlyingSampleThreshold              = None
        self.OutlierCorrelationPercentage          = None
        self.ImputationPackageName                 = None
        self.ImputationPackageType                 = None
        self.ImputationSettings                    = None
#       1c)
        self.ConcentrationTransformationFunction   = None
        self.PairwiseMetaboliteRatios              = None
#       1d)
        self.PositiveControl                       = None
        
#       Set values
        self.SetRefererenceMeanCVOfAllPlatesThreshold()
        self.SetMetaboliteMissingValueRateThreshold()
        self.SetMetaboliteOutlyingSampleThreshold()
        self.SetNOutlyingSampleThreshold()
        self.SetOutlierCorrelationPercentage()
        self.SetImputationPackageName()
        self.SetImputationPackageType()
        self.SetImputationSettings()
        self.SetConcentrationTransformationFunction()
        self.SetPairwiseMetaboliteRatios()
        self.SetPositiveControl()
        
    def SetRefererenceMeanCVOfAllPlatesThreshold(self):
        self.RefererenceMeanCVOfAllPlatesThreshold = 25.0 #%
        return
    def GetRefererenceMeanCVOfAllPlatesThreshold(self):
        return self.RefererenceMeanCVOfAllPlatesThreshold
    def SetMetaboliteMissingValueRateThreshold(self):
        self.MetaboliteMissingValueRateThreshold = 5.0 #%
        return
    def GetMetaboliteMissingValueRateThreshold(self):
        return self.MetaboliteMissingValueRateThreshold
    def SetMetaboliteOutlyingSampleThreshold(self):
        self.MetaboliteOutlyingSampleThreshold = 5.0 # is a factor with "unit" \sigma
        return
    def GetMetaboliteOutlyingSampleThreshold(self):
        return self.MetaboliteOutlyingSampleThreshold
    def SetNOutlyingSampleThreshold(self):
        self.NOutlyingSampleThreshold = 3 # Amount of independent samples (Sample Level)
        return
    def GetNOutlyingSampleThreshold(self):
        return self.NOutlyingSampleThreshold
    def SetOutlierCorrelationPercentage(self):
        self.OutlierCorrelationPercentage = 70.0 #%
        return
    def GetOutlierCorrelationPercentage(self):
        return self.OutlierCorrelationPercentage
    def SetImputationPackageName(self):
        self.ImputationPackageName = 'mice' # Multivariate imputation R-package 'mice'
        return
    def GetImputationPackageName(self):
        return self.ImputationPackageName
    def SetImputationPackageType(self):
        self.ImputationPackageType = 'R' # 'mice' is an R-package
        return
    def GetImputationPackageType(self):
        return self.ImputationPackageType
    def SetImputationSettings(self):
        self.ImputationSettings = None # AT THIS MOMTENT!!
        return
    def GetImputationSettings(self):
        return self.ImputationSettings
    def SetConcentrationTransformationFunction(self):
        self.ConcentrationTransformationFunction = scipy.log # scipy natural logarithm function
        return
    def GetConcentrationTransformationFunction(self):
        return self.ConcentrationTransformationFunction
    def SetPairwiseMetaboliteRatios(self):
        self.PairwiseMetaboliteRatios = 'a/b' # Note ln(a/b)=-ln(b/a) 
        return
    def GetPairwiseMetaboliteRatios(self):
        return self.PairwiseMetaboliteRatios
    def SetPositiveControl(self):
        self.PositiveControl  = 'SNP rs174547 (T/C) in FADS1'
        self.PositiveControl += 'for ration PC aa C36:3 / PC aa C36:4.\n'
        self.PositiveControl += 'Compare result with Illig et al. Nat. Genet, 42, 137-141, 2010.'
        return
    def GetPositiveControl(self):
        return self.PositiveControl