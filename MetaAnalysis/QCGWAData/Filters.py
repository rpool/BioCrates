
class Filters:
    def __init__(self):
        self.InitialArrayLength    = 0
        self.NDeletedDuplicateSNPs = 0
        self.MaxNDuplicateSNPs     = 0
        self.boDuplicateSNPWarning = False
        return

    def SetInitialArrayLength(self,
                              Value=int):
        self.InitialArrayLength = Value
        return

    def GetInitialArrayLength(self):
        return self.InitialArrayLength

    def SetNDeletedDuplicateSNPs(self,
                                 Value=int):
        self.NDeletedDuplicateSNPs = Value
        return

    def GetNDeletedDuplicateSNPs(self):
        return self.NDeletedDuplicateSNPs

    def RemoveDuplicateSNPs(self,
                            DCsDict={}):
        NRemoved = 0
        for Key in DCsDict.iterkeys():
            DuplicateIndexDict = DCsDict[Key].DataContainers['SNPID'].GetDuplicateIndexDict()
            for DCKey in DCsDict[Key].DataContainers.iterkeys():
                NRemoved = DCsDict[Key].DataContainers[DCKey].RemoveDuplicates(DuplicateIndexDict)

        self.SetNDeletedDuplicateSNPs(NRemoved)
        return DCsDict

    def SetMaxNDuplicateSNPs(self,
                             Value=int):
        self.MaxNDuplicateSNPs = Value
        return

    def GetMaxNDuplicateSNPs(self):
        return self.MaxNDuplicateSNPs

    def SetboDuplicateSNPWarning(self):
        if(self.GetMaxNDuplicateSNPs()>100):
            self.boDuplicateSNPWarning = False
        return

    def GetboDuplicateSNPWarning(self):
        return self.boDuplicateSNPWarning