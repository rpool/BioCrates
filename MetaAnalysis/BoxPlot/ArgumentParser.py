import argparse
import os

import Logger

#===============================================================================
# This module parses the command line options and generates a log
#===============================================================================

def LogArguments(Log=Logger,
                 ArgParser=argparse.ArgumentParser,
                 Arguments=argparse.Namespace):
#   Argument logging module
    ArgParser.print_help()
    ArgParser.print_help(Log.GetFileHandle())

#   Calculate max. of keylength for formatting
    MaxLen = 0
    for Key in vars(Arguments).iterkeys():
        MaxLen = max(MaxLen,len(Key))
    LogString =  '\n****\n'+\
                 'Used arguments:\n'+\
                 '---------------'
    print LogString
    Log.Write(LogString+'\n')

    FormatString = '{0:<'+str(MaxLen)+'}'
    LogString    = ''
    for Key,Value in vars(Arguments).iteritems():
        LogString += FormatString.format(Key)+': '+str(Value)+'\n'
    LogString += '****\n'
    print LogString
    Log.Write(LogString+'\n')
    return

def ParseArguments(Log=None):
#   Argument parsing module
    ArgumentParser = argparse.ArgumentParser(description=\
                                             'This python module can be used for box plotting GWA output values in a cross-cohort QC pipeline.')
    ArgumentParser.add_argument('-P',
                                '--protocolfile',
                                dest='ProtocolFile',
                                help='NAME: Name of file contains the necessary input and protocol values (input, format XML)',
                                metavar='NAME',
                                default=os.path.join(os.getcwd(),'Protocol.xml'))

    Arguments = ArgumentParser.parse_args()

    return ArgumentParser,\
           Arguments