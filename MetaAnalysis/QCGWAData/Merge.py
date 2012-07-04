import lxml.etree

def Merge(XmlObj=lxml.etree._ElementTree,
          SourceDCsDict={},
          DestDCsDict={}):

    for Column in XmlObj.getroot.find('MtbGWAColumns'):
        if(Column.find('FromExtraInfo')!=None):
            Tag = Column.tag
            Entry2IndexDict = SourceDCsDictGetEntry2IndexDict()
            for DKey in DestDCsDict.iterkeys():
                DC = DestDCsDict[DKey].DataContainers[Tag]

    for SKey in SourceDCsDict.iterkeys():


    return OutDCsDict