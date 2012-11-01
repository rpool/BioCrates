#! /usr/bin/env python

# python modules
import os
import fnmatch
import sys
import scipy
import re

# Homebrew modules
import Logger
import ArgumentParser
import File
import Plotting

def main(ExecutableName=str):

    #===========================================================================
    # START Initialization
    #===========================================================================

    ArgParser,\
    Arguments   = ArgumentParser.ParseArguments()

    Ext = ''
    Log = Logger.Logger(ExecutableName,
                        Ext)
    LogString  = '## START TIMESTAMP\n'
    LogString += str(Log.GetStartDate())+'\n'
    LogString += '## END TIMESTAMP'
    print LogString
    Log.Write(LogString+'\n')
    LogString = Log.GetStartLogString()
    print LogString
    Log.Write(LogString+'\n')
    ArgumentParser.LogArguments(Log,
                                ArgParser,
                                Arguments)
    #===========================================================================
    # END Initialization
    #===========================================================================

    #===========================================================================
    # START Do the work!
    #===========================================================================

    # Parse Arguments.ProtocolFile
    LogString = '++ Parsing \"'+Arguments.ProtocolFile+'\" ...'
    print LogString
    Log.Write(LogString+'\n')
    ProtocolFile = File.File(Name=Arguments.ProtocolFile,
                             boHeader=False)
    XmlProtocol  = ProtocolFile.ParseAsXml()
    LogString = '-- Done ...'
    print LogString
    Log.Write(LogString+'\n')

    ExcludedMtbs = []
    for Tag in XmlProtocol.getroot().find('ExcludedMtbs'):
        ExcludedMtbs.append(Tag.text)

    PhenotypePath       = XmlProtocol.getroot().find('PhenotypeInputPath').text.strip()
    PhenotypeExtension  = XmlProtocol.getroot().find('PhenotypeFileExt').text.strip()
    PhenotypeInputFiles = fnmatch.filter(os.listdir(PhenotypePath),'*.'+PhenotypeExtension)
    PhenotypeInputFiles = fnmatch.filter(PhenotypeInputFiles,'*_CHR1_*')
    PhenotypeFileDict   = {}
    for F in PhenotypeInputFiles:
        P = int(re.sub('PHE','',F.split('_')[-1].split('.')[0]))
        PhenotypeFileDict[P] = os.path.join(PhenotypePath,F)
    PhenotypeNameDict = {}
    PhenotypeExclDict = {}
    NColumns          = None
    for Key in PhenotypeFileDict.iterkeys():
        fr       = open(PhenotypeFileDict[Key],'r')
        LSplit   = fr.readline().strip().split()
        NColumns = len(LSplit)
        Name     = LSplit[-1]
        fr.close()
        PhenotypeNameDict[Key] = Name
        if(Name in ExcludedMtbs):
            PhenotypeExclDict[Key] = True
        else:
            PhenotypeExclDict[Key] = False
    PhenotypeArrayDict = {}
    for Key in PhenotypeFileDict.iterkeys():
        Array = scipy.loadtxt(fname=PhenotypeFileDict[Key],
                              dtype=str,
                              skiprows=2,
                              usecols=[NColumns-1],
                              unpack=True)
        FilterArray             = (Array!='NA')
        PhenotypeArrayDict[Key] = scipy.compress(FilterArray,Array).astype(float)

    for Key in PhenotypeFileDict.iterkeys():
        Name = PhenotypeNameDict[Key]
        if(not PhenotypeExclDict[Key]):
            print Key, Name, len(PhenotypeArrayDict[Key])

#    for i in range(len(PhenotypeInputFiles)):
#        p = str(i+1)

    #===========================================================================
    # END Do the work!
    #===========================================================================

    #===========================================================================
    # START Finalization
    #===========================================================================
    LogString = '\n**** Done :-)'
    print LogString
    Log.Write(LogString+'\n')
    LogString = Log.GetEndLogString()
    print LogString
    Log.Write(LogString+'\n')
    Log.Close()
    #===========================================================================
    # START Finalization
    #===========================================================================

    return

if(__name__=='__main__'):
    ExecutableName = os.path.abspath(__file__).split('/')[-1]
    main(ExecutableName)