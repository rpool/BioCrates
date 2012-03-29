class ConcentrationUnit:
    """
    Class containing the concentration unit.
    """
    
    def __init__(self,
                 String=str):
        self.MicroString = u'\xb5M'
        self.MicroFloat  = 1.0e-6
        self.Unit        = None
        self.UnitFloat   = None
        
        self.SetUnit(String)
          
    def SetUnit(self,
                String=str):
        self.Unit = String.split('[')[1]
        self.Unit = self.Unit.split(']')[0]
        
        if(self.Unit==self.MicroString):
            self.UnitFloat = self.MicroFloat
        return
        
    def GetUnit(self):
        return self.Unit
    
    def GetUnitFloat(self):
        return self.UnitFloat