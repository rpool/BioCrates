import lxml.etree
import scipy

def Merge(XmlObj=lxml.etree._ElementTree,
          SourceDCsDict={},
          DestDCsDict={},
          SourceColumnTag=str):

    Entry2IndexDictDict = {}
    for Key in SourceDCsDict.iterkeys():
        Entry2IndexDictDict[Key] = SourceDCsDict[Key].DataContainers[SourceColumnTag].GetEntry2IndexDict()

    for Column in XmlObj.getroot().find('MtbGWAColumns'):
        if((Column.find('GetFromExtraInfoFile')!=None) and
           (eval(Column.find('GetFromExtraInfoFile').text))):
            Tag = Column.tag
            for DKey in DestDCsDict.iterkeys():
                DataArray = []
                for Entry in DestDCsDict[DKey].DataContainers[SourceColumnTag].GetDataArray():
                    IndexInSourceDict = None
                    for SKey in Entry2IndexDictDict.iterkeys():
                        if(Entry2IndexDictDict[SKey].has_key(Entry)):
                            IndexInSourceDict = Entry2IndexDictDict[SKey][Entry]
                            break
                    SEntry = SourceDCsDict[SKey].DataContainers[Tag].GetDataArray()[IndexInSourceDict]
                    DataArray.append(SEntry)
                DataArray = scipy.array(DataArray)
                DestDCsDict[DKey].DataContainers[Tag].ReplaceDataArray(DataArray)

    return DestDCsDict