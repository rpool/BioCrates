import xlrd
import re

import ConcentrationUnit
import DataContainer
import Logger

def ParseExcellBook(ExcelDataFileName=str,
                    ExcelSheetName=str,
                    boRemoveCustomMetabolites=bool,
                    Log=Logger):
    ExcelBook  = xlrd.open_workbook(ExcelDataFileName)
    LogString =  '** Using Excel Sheet with name \"'+ExcelSheetName+'\" ...'
    print LogString
    Log.Write(LogString+'\n')
    ExcelSheet = ExcelBook.sheet_by_name(ExcelSheetName)
#   Get the info from the first row (this should contain "Concentration []")
    for Column in range(ExcelSheet.ncols):
        CellType = ExcelSheet.cell_type(0,Column)
        if(CellType!=xlrd.XL_CELL_EMPTY): # skip empty cells
                ConcUnit = ConcentrationUnit.ConcentrationUnit(ExcelSheet.cell_value(0,Column))
#   Get the info from the second row (this should contain the column names)
    DataDict = {}
    for Column in range(ExcelSheet.ncols):
        CellType = ExcelSheet.cell_type(1,Column)
        if(CellType!=xlrd.XL_CELL_EMPTY): # skip empty cells
            CellValue = ExcelSheet.cell_value(1,Column)
            DataDict[Column] = DataContainer.DataContainer()
            DataDict[Column].SetDataName(CellValue)
#   Get the info from the trird row (this should contain the column names)
    Counter = 0
    for Column in range(ExcelSheet.ncols):
        CellType = ExcelSheet.cell_type(2,Column)
        if(CellType!=xlrd.XL_CELL_EMPTY): # skip empty cells
            CellValue = ExcelSheet.cell_value(2,Column)
            if(CellValue!='Class'):
                DataDict[Column].SetMetaboliteClass(CellValue)
                DataDict[Column].SetMetaboliteName(DataDict[Column].GetDataName())
                Counter += 1
                DataDict[Column].SetMetaboliteIndex(Counter)

#   Get the info from the other rows
    for Row in range(3,ExcelSheet.nrows):
        boPassedLODInfoCell = False
        for Column in range(ExcelSheet.ncols):
            CellType = ExcelSheet.cell_type(Row,Column)
            if(ExcelSheet.cell_type(Row,Column)!=xlrd.XL_CELL_EMPTY): # skip empty cells
                CellValue = ExcelSheet.cell_value(Row,Column)
                if((CellType==xlrd.XL_CELL_TEXT) and
                   (re.search('LOD',CellValue))):
                    DataDict[Column].AppendToLODInfoArray(CellValue)
                    boPassedLODInfoCell = True
                elif(boPassedLODInfoCell):
                    DataDict[Column].AppendToLODValueArray(float(CellValue))
                elif((CellType==xlrd.XL_CELL_TEXT) or
                     (CellType==xlrd.XL_CELL_NUMBER)):
                    DataDict[Column].AppendToDataArray(CellValue)
#   Reomve custom metabolites?
    if(boRemoveCustomMetabolites):
        DelList = []
        for Key, Value in DataDict.iteritems():
            if(re.search('Cust.',str(Value.GetMetaboliteClass()))):
                DelList.append(Key)
        for Key in DelList:
            del DataDict[Key]

    return DataDict