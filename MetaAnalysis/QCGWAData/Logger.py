# -*- coding: utf-8 -*-
"""
This module handles the logging to stdout and to file.

.. moduleauthor:: René Pool <r.pool@vu.nl>

"""
import sys
import os
import pwd
import platform
import datetime
import re

#===============================================================================
# This module contains the basic Logger class.
# Its members and member modules should speak for themselves, if not a
# comment is provided.
#===============================================================================

class Logger:
    """
    Class **Logger**

    :param FileName: Name of the log file
    :type FileName: str
    :param Extension: Extension of the log file
    :type Extension: str -- default='.log'
    """
    def __init__(self,
                 FileName=str,
                 Extension=''):
        self.LogFileName = None
        if(Extension!=''):
            self.LogFileName = FileName+'.'+Extension+'.log'
        else:
            self.LogFileName = FileName+'.log'
        self.FileHandle     = None
        self.StartDate      = None
        self.StartLogString = None
        self.EndLogString   = None

        self.InitFileHandle()
        self.StartLog()

        return

    def InitFileHandle(self):
        """
        Initializes the log file handle.
        """
        ListDir        = os.listdir(os.getcwd())
#       The following is done to prevent previous log files to be overwritten.
#       NOTE: Counting goes until 9, so the *.log9 file will be overwritten in
#       case it is present!
        LogCounterList = [-1]
        for File in ListDir:
            if(re.search(self.GetFileName(),File)):
                LogCounterList.append(int(File[-1]))
        MaxCounter        = max(LogCounterList)
        self.LogFileName += str(MaxCounter+1)
        self.FileHandle   = open(self.LogFileName,'w')
        return
    def GetFileHandle(self):
        """
        Returns the log file handle.
        """
        return self.FileHandle

    def GetFileName(self):
        """
        Returns the log file name.
        """
        return self.LogFileName

    def GetStartDate(self):
        """
        Returns the log start date/time.
        """
        return self.StartDate

    def StartLog(self):
        """
        Writes the log header to the log file.
        """
        StartDate      = datetime.datetime.now()
        self.StartDate = StartDate

        self.StartLogString  = '\n'
        self.StartLogString += '#############\n'
        self.StartLogString += '# START LOG #\n'
        self.StartLogString += '#############\n'
        self.StartLogString += '#\n'
        self.StartLogString += '# STARTED AT:\n'
        self.StartLogString += '# date (yyyy-mm-dd) : '+str(StartDate.date())+'\n'
        self.StartLogString += '# time (hh:mm:ss)   : '+str(StartDate.time())+'\n'
        self.StartLogString += '# platform          : '+' '.join(platform.uname())+'\n'
        self.StartLogString += '# user              : '+str(pwd.getpwuid(os.getuid())[0])+'\n'
        self.StartLogString += '# path              : '+os.getcwd()+'\n'
        self.StartLogString += '# python version    : '+str(sys.version_info)+'\n'
        self.StartLogString += '# running module    : '+sys.argv[0]+'\n'
        if(os.path.islink(sys.argv[0])):
            self.StartLogString += '#  -> linking to    : '+os.path.realpath(sys.argv[0])+'\n'
        if(os.path.isdir(os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])),'.svn'))):
            try:
                import pysvn
                Client = pysvn.Client()
                self.StartLogString += '# svn revision      : '+str(Client.info(os.path.dirname(os.path.realpath(sys.argv[0]))).revision.number)+'\n'
                del Client
            except ImportError:
                fr    = open(os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])),'.svn','entries'),'r')
                Lines = fr.readlines()
                fr.close()
                for l in range(len(Lines)):
                    Line = Lines[l]
                    if(Line[:3]=='dir'):
                        break
                Version              = Lines[l+1].strip()
                self.StartLogString += '# svn revision      : '+Version+'\n'

        self.StartLogString += '\n'

        return

    def GetStartLogString(self):
        """
        Returns the log header string.
        """
        return self.StartLogString

    def EndLog(self):
        """
        Writes the log footer to the log file.
        """
        EndDate = datetime.datetime.now()

        self.EndLogString  = '\n'
        self.EndLogString += '###########\n'
        self.EndLogString += '# END LOG #\n'
        self.EndLogString += '###########\n'
        self.EndLogString += '#\n'
        self.EndLogString += '# FINISHED AT:\n'
        self.EndLogString += '# date (yyyy-mm-dd) : '+str(EndDate.date())+'\n'
        self.EndLogString += '# time (hh:mm:ss)   : '+str(EndDate.time())+'\n'
        self.EndLogString += '# elapsed (s)       : '+str((EndDate-self.StartDate).seconds) # Report total time.

        return

    def GetEndLogString(self):
        """
        Returns the log footer string.
        """
        self.EndLog()
        return self.EndLogString

    def Write(self,
              String=str):
        """
        Writes a string to the log file.

        :param String: The string to be written to the log file
        :type String: str
        :returns: nothing
        """
        self.FileHandle.write(String)
        return

    def Close(self):
        """
        Closes the log file handle.
        """
        self.EndLog()
        self.FileHandle.close()
        return
