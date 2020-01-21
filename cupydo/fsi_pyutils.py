#! /usr/bin/env python
# -*- coding: utf8 -*-

''' 

Copyright 2018 University of Liège

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 

fsi_pyutils.py
Python utilities.
Authors: L. PAPELEUX

'''
from __future__ import print_function

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

from future import standard_library
standard_library.install_aliases()
import os, os.path, sys, platform
import tempfile

# ----------------------------------------------------------------------
#  Utilities
# ----------------------------------------------------------------------

def quit():
    sys.exit()

def fileToModule(file, verb=True):
    """
    convert a path/file to a module name (e.g. apps/qs/cont2.py => apps.qs.cont2)
    """
    file=os.path.abspath(file)
    if verb: print('file=',file)
    for dirname in sys.path:
        dirname = os.path.abspath(dirname)
        if verb: print('module in', dirname, '?')

        if isUnix(): 
            common = os.path.commonprefix( (file, dirname) )
        else:
            dirname = dirname.lower()
            common = os.path.commonprefix( (file.lower(), dirname) )
        
        if common == dirname:
            strip = file[len(dirname):]
            if strip[0]==os.path.sep:  # enl�ve le slash si pr�sent dans sys.path
                strip=strip[1:]
            strip = os.path.splitext(strip)[0]
            if strip.find('.')!=-1: # dir contains dots
                continue
            strip = strip.replace(os.path.sep,'.')
            # on verifie qu'il y a des __init__.py tout le long
            dirs = strip.split('.'); dirs.pop()
            dir=common
            for sub in dirs:
                dir = os.path.join(dir,sub)
                init_py = os.path.join(dir,"__init__.py")
                if not os.path.isfile(init_py):
                    #if verb: print "creating", init_py
                    touchFile = open(init_py, 'w')
                    touchFile.close()
            if verb: 
                print('YES')
                print('module=', strip)
            return strip
            break
        else:
            if verb: print('NO')
    return ''

def isUnix():
    uname = platform.uname()
    return not (uname[0] == 'Windows' or uname[2] == 'Windows')
    
def getWinTxtEditor():
    if os.path.isfile('c:\\\\Program Files\\\\JGSoft\\\\EditPadPro5\\\\EditPadPro.exe') :
        return 'c:\\\\Program Files\\\\JGSoft\\\\EditPadPro5\\\\EditPadPro.exe'
    elif os.path.isfile('c:\\\\Program Files\\\\JGSoft\\\\EditPadPro6\\\\EditPadPro.exe') :
        return 'c:\\\\Program Files\\\\JGSoft\\\\EditPadPro6\\\\EditPadPro.exe'
    elif os.path.isfile('c:\\\\Program Files (x86)\\\\JGSoft\\\\EditPadPro6\\\\EditPadPro.exe') :
        return 'c:\\\\Program Files (x86)\\\\JGSoft\\\\EditPadPro6\\\\EditPadPro.exe'
    else :
        return 'c:\\\\windows\\\\notepad.exe'
        
def hasEditPad():
    return os.path.isfile("c:\\\\Program Files\\\\JGSoft\\\\EditPadPro5\\\\EditPadPro.exe")

def pythonPathSep():
    if isUnix():
        return ':'
    else:
        return ';'

def cmdInPath(cmd, othdirs):
    """
    find an exec in the PATH or othdirs. Return the path to the exec
    """
    foundfile=None
    for p in othdirs + os.environ['PATH'].split(pythonPathSep()):
        file1 = os.path.join(p, cmd)
        file2 = file1+'.exe'
        for f in (file1,file2):
            if os.path.isfile(f):
                foundfile=f
                return foundfile
    return foundfile
    
def canCreateFile(folder_path):
    try:
        tempfile.TemporaryFile(dir=folder_path)
        return True
    except OSError:
        return False

def canCreateFolder(folder_path):
    try:
        name = tempfile.mkdtemp(dir=folder_path)
        os.rmdir(name)
        return True
    except OSError:
        return False

def printMem(indent=""):    
    if isUnix():
        import subprocess
        #
        procId = os.getpid()
        res = subprocess.getoutput('cat /proc/%s/status' % procId).split('\n')
        status = dict()
        for i in res:
            if i != '':
                res2 = i.split(":\t")                
                status[res2[0]]=res2[1]    
        print(indent+"VmSize: %s  VmRSS: %s  VmData: %s " % (status['VmSize'], status['VmRSS'], status['VmData']))
    else:
        try:
            import win32process
            import win32con
            import win32api
            procId = win32process.GetCurrentProcessId()
            han    = win32api.OpenProcess(win32con.PROCESS_QUERY_INFORMATION|win32con.PROCESS_VM_READ, 0, procId)
            procmeminfo  = win32process.GetProcessMemoryInfo(han)     
            workMem      = (procmeminfo["WorkingSetSize"]/1024.)   
            peakWork     = (procmeminfo["PeakWorkingSetSize"]/1024.)
            pageFile     = (procmeminfo["PagefileUsage"]/1024.)
            peakPageFile = (procmeminfo["PeakPagefileUsage"]/1024.)
            print(indent+"WorkMem: %sK  PeakMem: %sK  PageFile: %sK  PeakPageFile: %sK" % (workMem,peakWork,pageFile,peakPageFile))
        except:
            print("install pywin32 to be able to compute process memory")
    
