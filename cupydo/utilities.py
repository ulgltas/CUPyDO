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

utilities.py
Common utilities (MPI functions, timer, ...) for CUPyDO.
Authors : David THOMAS, Marco Lucio CERQUAGLIA, Romain BOMAN, Adrien CROVATO

'''

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

from math import *
import numpy as np
import scipy as sp
import os, os.path, sys
import time as tm
import json

np.set_printoptions(threshold=sys.maxsize)

# global vars (underscore prevent them to be imported with "from module import *")
_theModule  = None
_theWDir    = None # workspace directory
_theWDirRoot = os.getcwd()  # base directory du workspace

# ----------------------------------------------------------------------
#  Utilities
# ----------------------------------------------------------------------

def solve_upper_triangular_mod(U, y, toll):

    # 'Modified' backward solve Ux = y with U upper triangular. 'Modified' because if the value of a diagonal element is close to zero (i.e. <= toll) the corresponding element in the solution is set to zero

    n = y.shape[0]
    x = np.zeros((n,))

    for i in range (n - 1, -1, -1):
        if abs(U[i, i]) <= toll:
            x[i] = 0.
        else:
            x[i] = y[i]
            for j in range (i + 1, n):
                  x[i] -= U[i, j] * x[j]
            x[i] /= U[i, i]
    return x.T

def QRfiltering(V, W, toll):
    
    while True:
        n = V.shape[0]
        s = V.shape[1]
        
        flag = True
        
        Q, R = sp.linalg.qr(V, mode='economic')
        
        for i in range(0, s): # REMEMBER: range(a, b) starts at a but ends at b-1!
            if abs(R[i,i]) < toll*np.linalg.norm(R, 2):
                V = np.delete(V, i, 1)
                W = np.delete(W, i, 1)
                if i == s-1:
                    flag = False
                else:
                    flag = True
                break
        if i >= s-1 and flag == True:
            return (Q, R, V, W)

def QRfiltering_mod(V, W, toll):
    
    while True:
        n = V.shape[0]
        s = V.shape[1]
        Q = np.zeros((n, s))
        R = np.zeros((s, s))
        
        flag = True
        
        i = 0
        V0 = V[:,i]
        R[i,i] = np.linalg.norm(V0, 2)
        Q[:,i] = np.dot(V0, 1.0/R[i,i])
        
        for i in range(1, s): # REMEMBER: range(a, b) starts at a but ends at b-1!
            vbar = V[:,i]
            for j in range(0, i): # REMEMBER: range(a, b) starts at a but ends at b-1!
                R[j,i] = np.dot(Q[:,j].T,vbar)
                vbar = vbar - np.dot(Q[:,j], R[j,i])
            if np.linalg.norm(vbar, 2) < toll*np.linalg.norm(V[:,i], 2):
                V = np.delete(V, i, 1)
                W = np.delete(W, i, 1)
                if i == s-1:
                    flag = False
                else:
                    flag = True
                break
            else:
                R[i,i] = np.linalg.norm(vbar, 2)
                Q[:,i] = np.dot(vbar, 1.0/R[i,i])
        if i >= s-1 and flag == True:
            return (Q, R, V, W)

# ------------------------------------------------------------------------------

def parseArgs():
    """Parse arguments specific to CUPyDO
    """
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('file', nargs=1, help='Python input file')
    parser.add_argument('-n', dest='n', help='Number of threads', default=1, type=int)
    parser.add_argument('--nogui', action='store_true', help='Disable GUI', default=False)
    parser.add_argument('-k', nargs=1, help='Number of threads')
    return parser.parse_args()

def setDirs(fpath):
    """Add file dir to python path and setup the working directory
    """
    import os, sys
    # get MPI status
    haveMpi, comm, rank, size = getMpi()
    # create workspace path
    sys.path.append(os.path.dirname(fpath)) # [RB] !this folder is can be a subfolder of one of the folders in the pythonpath
    base = os.path.abspath(os.path.join(os.path.dirname(os.path.dirname(__file__)), ' ')) # find base directory
    mpiPrint(base, comm)
    common = os.path.commonprefix((fpath, base)) # find common part of testname ad base name
    resdir = os.path.splitext(fpath[len(common):].replace(os.sep,"_"))[0] # common part, change seprator to underscore and remove ".py"
    wdir=os.path.join('workspace', resdir) # that is our workspace!
    # create workspace (master process)
    if rank == 0:
        if not os.path.isdir(wdir):
            print("creating", wdir)
            os.makedirs(wdir)

        if size > 1: # send sync to slaves
                for i in range(1, size):
                    comm.send(1, dest=i, tag=11)

        if comm:
            comm.barrier()
    # block slave processes
    else:
        comm.barrier()
    # change to workspace
    mpiPrint('changing to ' + wdir, comm)
    os.chdir(wdir)

def getMpi():
    """
    Get MPI parameters
    """
    from ccupydo import CMpi
    cmpi = CMpi()
    if cmpi.haveMPI:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        myid = comm.Get_rank()
        numberPart = comm.Get_size()
    else:
        comm = None
        myid = 0
        numberPart = 1
    return cmpi.haveMPI, comm, myid, numberPart

# ----------------------------------------------------------------------
#    MPI Functions
# ----------------------------------------------------------------------

def mpiPrint(message, mpiComm = None):
    """
    Description.
    """

    if mpiComm != None:
        myid = mpiComm.Get_rank()
    else:
        myid = 0

    if myid == 0:
        print(message)

    mpiBarrier(mpiComm)

def mpiBarrier(mpiComm = None):
    """
    Description.
    """

    if mpiComm != None:
        mpiComm.barrier()

def mpiAllReduce(mpiComm = None, value = 0):
    """
    Description.
    """
    sendBuff = np.array(value)
    if mpiComm != None:
        from mpi4py import MPI
        rcvBuff = np.zeros(1, dtype=type(value))
        mpiComm.Allreduce(sendBuff, rcvBuff, MPI.SUM)
        return rcvBuff[0]
    else:
        return value

def mpiAllGather(mpiComm = None, value = 0):
    """
    Description
    """

    sendBuff = np.array(value)
    if mpiComm != None:
        mpiSize = mpiComm.Get_size()
        rcvBuff = np.zeros(mpiSize, dtype=type(value))
        mpiComm.Allgather(sendBuff, rcvBuff)
        return rcvBuff
    else:
        return sendBuff

def mpiGatherv(sendBuff, localSize, globalSize, mpiComm = None, rootProcess=0):
    """
    Des.
    """

    #sendBuff must be a numpy array

    rcvBuff = None

    if mpiComm != None:
        from mpi4py import MPI
        myid = mpiComm.Get_rank()
        mpiSize = mpiComm.Get_size()
        if myid == rootProcess:
            rcvBuff = np.zeros(globalSize)
        sendBuffNumber = np.array([localSize], dtype=int)
        rcvBuffNumber =np.zeros(mpiSize, dtype=int)
        mpiComm.Allgather(sendBuffNumber, rcvBuffNumber)
        counts = tuple(rcvBuffNumber)
        displ = np.zeros(mpiSize, dtype=int)
        for ii in range(rcvBuffNumber.shape[0]):
            displ[ii] = rcvBuffNumber[0:ii].sum()
        displ = tuple(displ)
        mpiComm.Gatherv(sendBuff, [rcvBuff, counts, displ, MPI.DOUBLE], root=rootProcess)
        return rcvBuff
    else:
        return sendBuff

def mpiGatherInterfaceData(interfData, globalSize, mpiComm = None, rootProcess = 0):
    """
    Des.
    """
    interfData_Gat = []

    if mpiComm != None :
        localSize = interfData.getDataArray(0).shape[0]
        for iDim in range(interfData.nDim):
            array_Gat = mpiGatherv(interfData.getDataArray(iDim), localSize, globalSize, mpiComm, 0)
            interfData_Gat.append(array_Gat)
    else:
        for iDim in range(interfData.nDim):
            interfData_Gat.append(interfData.getData(iDim))

    return interfData_Gat

# ----------------------------------------------------------------------
#   Timer class
# ----------------------------------------------------------------------

class Timer(object):
    """
    Description
    """

    def __init__(self):
        """
        Des.
        """

        self.startTime = 0.0
        self.stopTime = 0.0
        self.elapsedTime = 0.0
        self.cumulTime = 0.0
        self.isRunning = False

    def start(self):
        """
        Des.
        """

        if not self.isRunning:
            self.startTime = tm.time()
            self.isRunning = True
        else:
            print('[WARNING] Calling Timer.start() but the Timer is already running !')

    def stop(self):
        """
        Des.
        """

        if self.isRunning:
            self.stopTime = tm.time()
            self.elapsedTime = self.stopTime - self.startTime
            self.isRunning = False

    def cumul(self):
        """
        Des.
        """

        self.cumulTime += self.elapsedTime
        self.elapsedTime = 0.0

    def getElapsedTime(self):
        """
        Des.
        """

        return self.elapsedTime

    def getCumulTime(self):
        """
        Des.
        """

        return self.cumulTime

# ----------------------------------------------------------------------
#   Solver path class
# ----------------------------------------------------------------------

class solverPath(object):
    def __init__(self):

        # Path to the base folder of CUPyDO

        baseFolder = __file__
        for i in range(3): baseFolder = os.path.split(baseFolder)[0]
        self.addedPaths = {}

        # Look for cupydoExtProgs.json in CUPyDO folder

        try:
            with open(os.path.join(baseFolder,'.cupydoExtProgs.json')) as file:
                self.solverPaths = json.load(file)

        # LIf not found, looks for cupydoExtProgs.json in home folder
        
        except:
            try:
                with open(os.path.join(os.path.expanduser('~'),'.cupydoExtProgs.json')) as file:
                    self.solverPaths = json.load(file)
            except: self.solverPaths = {}

        # Relative paths for solvers located next to CUPyDO

        self.basePaths = {}
        self.basePaths['Metafor'] = ['Metafor/oo_metaB/bin','Metafor/oo_meta','Metafor/linuxbin']
        self.basePaths['RBMI'] = ['NativeSolid/bin']
        self.basePaths['Modal'] = ['modali']
        self.basePaths['Flow'] = ['waves']
        self.basePaths['Pfem'] = ['PFEM','waves']
        self.basePaths['SU2'] = ['SU2/build/bin']
        self.basePaths['VLM'] = ['VLM']
        self.basePaths['Pfem3D'] = ['PFEM3D/build/bin']
        self.basePaths['pyBeam'] = ['pyBeam/build/bin']
        
        # Makes the absolute paths for solvers next to CUPyDO

        for key,value in self.basePaths.items():
            for i in range(len(value)):

                value[i] = os.path.join(baseFolder,value[i])
                release = os.path.join(value[i],'Release')
                if os.path.isdir(release): value[i] = release

    # Adds the path to the solver in sys.path

    def add(self, solverName):

        # Check next to the solver first, then in the Json

        try: pathList = self.solverPaths[solverName]
        except:
            try: pathList = self.basePaths[solverName]
            except: raise Exception('%s is not found' % solverName)
        
        if not isinstance(pathList,list): pathList = [pathList]

        # Memorize the paths that have been added under that name

        if solverName in self.addedPaths:
            self.addedPaths[solverName] += pathList
        else: self.addedPaths[solverName] = pathList

        # Actually adds the paths to the system

        for path in pathList:
            if os.path.isdir(path):
                sys.path.append(path)
                print('Add to path :',path)
            else:
                raise Exception('%s not found' % path)

    # Remove the path to the solver from sys.path

    def remove(self, solverName):
        
        try: pathList = self.addedPaths[solverName]
        except: raise Exception('%s has not been added' % solverName)

        if not isinstance(pathList,list): pathList = [pathList]

        for path in pathList:
            sys.path.remove(path)
            print('Remove from path :',path)