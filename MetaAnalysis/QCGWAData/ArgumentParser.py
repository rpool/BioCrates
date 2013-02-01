# -*- coding: utf-8 -*-
"""
This module handles the argument parsing.

.. moduleauthor:: René Pool <r.pool@vu.nl>

"""

import argparse
import os

import Logger

#===============================================================================
# This module parses the command line options and generates a log
#===============================================================================

def LogArguments(Log=Logger,
                 ArgParser=argparse.ArgumentParser,
                 Arguments=argparse.Namespace):
    """
    This function handles the argument logging to stdout and to file.

    :param Log: A Logger instance
    :type Log: :class:Logger
    :param ArgParser: An argparse.ArgumentParser instance
    :type ArgParser: argparse.ArgumentParser
    :param Arguments: Arguments read from stdin
    :type Arguments: argparse.Namespace
    :returns: nothing
    """
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
    """
    This function handles the argument logging to stdout and to file.

    :param Log: A Logger instance (dummy argument)
    :type Log: None
    :param ArgParser: An argparse.ArgumentParser instance
    :type ArgParser: argparse.ArgumentParser
    :param Arguments: Arguments read from stdin
    :type Arguments: argparse.Namespace
    :returns: **ArgumentParser** (an argparse.ArgumentParser instance), **Arguments** (an argparse.Namespace instance)
    """
#   Argument parsing module
    ArgumentParser = argparse.ArgumentParser(description=\
                                             'This python module can be used for the initial QC of GWA output files.')
    ArgumentParser.add_argument('-P',
                                '--protocolfile',
                                dest='ProtocolFile',
                                help='NAME: Name of file contains the necessary input and protocol values (input, format XML)',
                                metavar='NAME',
                                default=os.path.join(os.getcwd(),'Protocol.xml'))

    Arguments = ArgumentParser.parse_args()

    return ArgumentParser,\
           Arguments