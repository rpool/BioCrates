class BioCratesAnalyticalRanges:
    def __init__(self,
                 FileName=str):
        fr = open(FileName,'r')
        FMem = fr.readlines()
        fr.close()
        self.LODValueDict  = {}
        self.LLOQValueDict = {}
        self.ULOQValueDict = {}
        for i in range(1,len(FMem)):
            Line      = FMem[i]
            IndexList = []
            IndexList.append(0)
            for j in range(len(Line)):
                Char = Line[j]
                if(Char==','):
                    IndexList.append(j)
            IndexList.append(len(Line))
            Name = None
            for i in range(1,len(IndexList)):
                StartIndex = IndexList[i-1]
                if(i==1):
                    Name = Line[StartIndex:IndexList[i]].strip()
                elif(i==4):
                    LODValue = float(Line[StartIndex+1:IndexList[i]].strip())
                    self.LODValueDict[Name] = LODValue
                elif(i==5):
                    LLOQValue = Line[StartIndex+1:IndexList[i]].strip()
                    if(len(LLOQValue)==0):
                        LLOQValue = None
                    else:
                        LLOQValue = float(LLOQValue)
                    self.LLOQValueDict[Name] = LLOQValue
                elif(i==6):
                    ULOQValue = Line[StartIndex+1:IndexList[i]].strip()
                    if(len(ULOQValue)==0):
                        ULOQValue = None
                    else:
                        ULOQValue = float(ULOQValue)
                    self.ULOQValueDict[Name] = ULOQValue
    def GetLODValueDict(self):
        return self.LODValueDict
    def GetLLOQValueDict(self):
        return self.LLOQValueDict
    def GetULOQValueDict(self):
        return self.ULOQValueDict