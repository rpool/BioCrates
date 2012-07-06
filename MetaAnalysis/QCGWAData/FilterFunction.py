import operator
import scipy

class FilterFunction:
    def __init__(self,
                 OperatorString=str,
                 CompareString=str,
                 CompareType=str):
        self.OperatorDict   = None
        self.TypeDict       = None
        self.OperatorString = OperatorString
        self.CompareString  = CompareString
        self.CompareType    = CompareType

        self.SetOperatorDict()
        self.SetTypeDict()
        return

    def Run(self,
            DataArray=scipy.array):
        Operator = self.GetOperatorDict()[self.OperatorString]
        Type     = self.GetTypeDict()[self.CompareType]

        FilterArray = Operator(DataArray.astype(Type),Type(self.CompareString))

        return FilterArray

    def SetTypeDict(self):
        self.TypeDict          = {}
        self.TypeDict['float'] = float
        self.TypeDict['int']   = int
        self.TypeDict['str']   = str
        return

    def GetTypeDict(self):
        return self.TypeDict

    def SetOperatorDict(self):
        self.OperatorDict       = {}
        self.OperatorDict['NE'] = operator.ne
        self.OperatorDict['EQ'] = operator.eq
        self.OperatorDict['GT'] = operator.gt
        self.OperatorDict['LT'] = operator.lt
        self.OperatorDict['GE'] = operator.ge
        self.OperatorDict['LE'] = operator.le
        return

    def GetOperatorDict(self):
        return self.OperatorDict
