import re
import gzip
import scipy
import sys

import Logger

#===============================================================================
# This module contains the basic DataContainer and DataContainers classes.
# Their members and member modules should speak for themselves, if not a
# comment is provided.
#===============================================================================

class DataContainer:
    def __init__(self):
        self.DataName  = None
        self.DataArray = None
        return

    def SetDataName(self,
                    Name=str):
        self.DataName = Name
        return

    def GetDataName(self):
        return self.DataName

    def InitDataArray(self):
        self.DataArray = []
        return

    def AppendToArray(self,
                      Entry=str):
        self.DataArray.append(Entry)
        return

    def CastDataArrayToScipy(self):
        self.DataArray = scipy.array(self.GetDataArray())

    def GetDataArray(self):
        return self.DataArray

class DataContainers:
    def __init__(self):
        self.DataContainers = {}
        self.Names2Columns  = {}
        self.Columns2Names  = {}
        return