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
                                             'This python module can be used for generating QQ plots in p-value and/or score mode.')
    ArgumentParser.add_argument('-M',
                                '--mtbnamefile',
                                dest='MtbNameFile',
                                help='PATH: Name of file that lists the GWA output files to be analyzed (input)',
                                metavar='PATH',
                                default=os.path.join(os.getcwd(),'MtbNames.txt'))
    ArgumentParser.add_argument('-p',
                                '--GWAdatapath',
                                dest='GWADataPath',
                                help='PATH: Name of path that contains the GWA output files to be analyzed (.npy input)',
                                metavar='PATH',
                                default=os.path.join(os.getcwd(),'PyBinData'))
    ArgumentParser.add_argument('-m',
                                '--qqmodes',
                                dest='QQModes',
                                help='PROPERTY: Determines the QQModes you want to run (input): '+\
                                     '\"P\"  sets the \"p-value\" mode; '+\
                                     '\"S\"  sets the \"score\" mode; '+\
                                     '\"PS\" sets them both; '+\
                                     '\"None\" uses no strata.',
                                metavar='PROPERTY',
                                default='P',
                                choices=['P','S','PS','None'])
    ArgumentParser.add_argument('-r',
                                '--PDFReport',
                                dest='boGeneratePdfReport',
                                help='FLAG: Generate a QQ analysis report in pdf format (input).',
                                action='store_true',
                                default=False)
    ArgumentParser.add_argument('-i',
                                '--ISqLt75',
                                dest='boFilterOnIQsLt75',
                                help='FLAG: Filter data on heterogeneity I^2 value .LT. 75.0%% (input).',
                                action='store_true',
                                default=False)

    Arguments = ArgumentParser.parse_args()

    return ArgumentParser,\
           Arguments