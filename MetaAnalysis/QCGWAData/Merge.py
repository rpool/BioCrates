import lxml.etree
import scipy
import copy

def MergeWithExtraInfo(XmlObj=lxml.etree._ElementTree,
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
                DataArray = copy.deepcopy(DestDCsDict[DKey].DataContainers[Tag].GetDataArray())
                DataArray = DataArray.tolist()
                for i in range(len(DestDCsDict[DKey].DataContainers[SourceColumnTag].GetDataArray())):
                    Entry = DestDCsDict[DKey].DataContainers[SourceColumnTag].GetDataArray()[i]
                    IndexInSourceDict = None
                    for SKey in Entry2IndexDictDict.iterkeys():
                        if(Entry2IndexDictDict[SKey].has_key(Entry)):
                            IndexInSourceDict = Entry2IndexDictDict[SKey][Entry]
                            break
                    if(IndexInSourceDict!=None):
                        SEntry       = SourceDCsDict[SKey].DataContainers[Tag].GetDataArray()[IndexInSourceDict]
                        DataArray[i] = SEntry
                DataArray = scipy.array(DataArray)
                DestDCsDict[DKey].DataContainers[Tag].ReplaceDataArray(DataArray)

    return DestDCsDict