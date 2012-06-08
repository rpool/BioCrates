import argparse
import os

import Logger

def LogArguments(Log=Logger,
                 ArgParser=argparse.ArgumentParser,
                 Arguments=argparse.Namespace):
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

    FormatString = '{:<'+str(MaxLen)+'}'
    LogString    = ''
    for Key,Value in vars(Arguments).iteritems():
        LogString += FormatString.format(Key)+': '+str(Value)+'\n'
    LogString += '****\n'
    print LogString
    Log.Write(LogString+'\n')
    return

def ParseArguments(Log=None):
    ArgumentParser = argparse.ArgumentParser(description=\
                                             'This python module can be used for generating QQ plots in p-value and/or score mode.')
    ArgumentParser.add_argument('-f',
                                '--gwafiles',
                                dest='GwaFiles',
                                help='PATH: Name of file that lists the GWA output files to be analyzes (input)',
                                metavar='PATH',
                                default=os.path.join(os.getcwd(),'GWAFiles.txt'))
    ArgumentParser.add_argument('-m',
                                '--qqmodes',
                                dest='QQModes',
                                help='PROPERTY: Determines the QQModes you want to run (input): '+\
                                     '\"P\"  sets the \"p-value\" mode; '+\
                                     '\"S\"  sets the \"score\" mode;  '+\
                                     '\"I\"  sets the \"score\" mode; '+\
                                     '\"XY\", where X=["P","S"] and Y=["S","I"] and X!=Y, sets two modes;'
                                     '\"PSI\" sets them all.',
                                metavar='PROPERTY',
                                default='P',
                                choices=['P','S','PS','PI','SI','PSI'])
    ArgumentParser.add_argument('-x',
                                '--xproperty',
                                dest='XProperty',
                                help='STRING: Property to be displayed on the x-axis (input)',
                                metavar='STRING',
                                default='pos')
    ArgumentParser.add_argument('-y',
                                '--yproperty',
                                dest='YProperty',
                                help='STRING: Property to be displayed on the y-axis (input)',
                                metavar='STRING',
                                default='PHE')

    Arguments = ArgumentParser.parse_args()

    return ArgumentParser,\
           Arguments