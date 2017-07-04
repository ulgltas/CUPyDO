#!/usr/bin/env python

## \file FSICoupler.py
#    \brief Description.
#    \author 
#    \version BETA
#
# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

from math import *
import numpy as np
import scipy as sp
from scipy import spatial
import scipy.sparse.linalg as splinalg
import os, os.path, sys, string
import time as tm

import socket, fnmatch
import fsi_pyutils

import copy

np.set_printoptions(threshold=np.nan)

# global vars (underscore prevent them to be imported with "from module import *")
_theModule  = None
_theWDir    = None # workspace directory
_theWDirRoot = os.getcwd()  # base directory du workspace

def solve_upper_triangular_mod(U, y, toll):
    
    # 'Modified' backward solve Ux = y with U upper triangular. 'Modified' because if the value of a diagonal element is close to zero (i.e. <= toll) the corresponding element in the solution is set to zero 
    
    n = y.shape[0]
    x = np.zeros((n,1))
    
    for i in range (n - 1, -1, -1):
        if fabs(U[i, i]) <= toll:
            x[i] = 0.
        else:
            x[i] = y[i];
            for j in range (i + 1, n):
                  x[i] -= U[i, j] * x[j]
            x[i] /= U[i, i]
    return x

# ------------------------------------------------------------------------------

def load(fsiTxt, mpi_opt, com, my_id, number_part):
    """load a module and make it the current one
    """
    global _theModule
    
    if not isinstance(fsiTxt, str): raise Exception("argument must be a string!")
    if _theModule: raise Exception("no more than one module may be loaded!")

    if fsiTxt=="__main__": # on est dans un script => getpfem est/sera dans __main__
        _theModule = sys.modules['__main__']
    else: 
        if os.path.isfile(fsiTxt): # c'est une nom de fichier => on convertit en nom de module
            module = pyutils.fileToModule(fsiTxt, verb=False)
            if not module: raise Exception('file "%s" is not a reachable python module!' % fsiTxt)
            fsiTxt = module
        # ici, on a un nom de module...
        exec("import %s" % fsiTxt) # peut on le faire plus tard?
        exec("globals()['_theModule'] = %s" % fsiTxt)
        print "module '%s' loaded!" % fsiTxt
        #print '_theModule', _theModule
        
    setTheWDir(fsiTxt)
    setDir(_theWDir)
    _chkWorkspace()
    
    if mpi_opt == True:
        from mpi4py import MPI  # MPI is initialized from now by python and can be continued in C++ !
        com = MPI.COMM_WORLD
        my_id = com.Get_rank()
        number_part = com.Get_size()
    else:
        com = None
        my_id = 0
        number_part = 1
    
    return _theModule

def setTheWDir(fsiTxt): 
    global _theWDir
    _theWDir = defaultWDir(fsiTxt)
    # on fait un chdir si ce rep existe!
    if os.path.isdir(_theWDir):
        setDir(_theWDir)
        
    return _theWDir

def defaultWDir(moduleTxt):
    global _theWDirRoot
    return os.path.join(_theWDirRoot, os.path.join("workspace", moduleTxt.replace('.','_') )) # default WDir for the module

def setDir(wdir):  # (avant instance) - change le rep courant
    global _theWDirRoot
    """ change the current working directory
    """
    if wdir:
        if not os.path.isabs(wdir):
            os.path.join(_theWDirRoot, wdir)
        try:
            if not os.path.isdir(wdir):
                os.makedirs(wdir)
            os.chdir(wdir)
        except OSError, e:
            print "OSError : ", e
            # check UAC (Vista) or root-protected dirs (Unix)    
            text = 'Directory "%s" may not be created.\n' % wdir
            text += 'Disable UAC (vista) or ask for root privilege (Unix) if you still want to work here.\n'
            text += 'Otherwise, rebase to another directory!' 
            try: wrap.GUILink().warningMsg(text)
            except: pass
            raise Exception('directory "%s" cannot be created' % wdir)    
    global _theWDir
    _theWDir = os.getcwd() # stoque un chemin absolu meme si ca a merde.
    print "_theWDir = ", _theWDir

def _chkWorkspace():
    """ Check/Clean workspace
    """
    flist=[]
    for file in os.listdir('.'):
        if os.path.isdir(file):
            raise Exception("Your workspace contains one or more directories!")
        elif os.path.isfile(file):
            for sk in ['*.fdb','*.conf','*.msh']: # files to keep
                if fnmatch.fnmatch(os.path.basename(file), sk):
                    break
            else:
                flist.append(os.path.abspath(file))
    
    if flist:
        answer = True
        try: answer = wrap.GUILink().askDelete(flist)
        except: pass
        if answer:
            for file in flist:
                try:
                    os.remove(file)
                except:
                    print "Unable to remove", file
        else:
            raise Exception("meta() cancelled!")


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
        return sendBuff

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
#   FlexInterfaceData class
# ----------------------------------------------------------------------

class FlexInterfaceData:
    """
    Description
    """

    def __init__(self, val_nPoint, val_nDim, mpiComm=None):
        """
        Des.
        """

        self.mpiComm = mpiComm
        self.nPoint = val_nPoint
        self.nDim = val_nDim

        self.dataContainer = []

        if mpiComm != None:
            for iDim in range(self.nDim):
                from petsc4py import PETSc
                data = PETSc.Vec().create(self.mpiComm)
                data.setType('mpi')
                data.setSizes(self.nPoint)
                data.set(0.0)
                self.dataContainer.append(data)
            self.myid = self.mpiComm.Get_rank()
            self.mpiSize = self.mpiComm.Get_size()
            startIndex , stopIndex = self.dataContainer[0].getOwnershipRange()
            self.indexing = self.mpiComm.allgather((startIndex , stopIndex))
        else:
            for iDim in range(self.nDim):
                data = np.zeros(self.nPoint, dtype=float)
                self.dataContainer.append(data)
            self.myid = 0
            self.mpiSize = 1

    def getnPoint(self):
        """
        Des.
        """

        return self.nPoint

    def getDim(self):
        """
        Des.
        """

        return self.nDim

    def setValues(self, iDim, indices_list, values_array):
        """
        Des.
        """

        if self.mpiComm != None:
            self.dataContainer[iDim].setValues(indices_list, values_array)
        else:
            self.dataContainer[iDim][indices_list] = values_array

    def setAllValues(self, iDim,value):
        """
        Des.
        """

        if self.mpiComm != None:
            self.dataContainer[iDim].set(value)
        else:
            self.dataContainer[iDim].fill(value)

    def getDataContainer(self):
        """
        Des.
        """

        return self.dataContainer

    def getData(self, iDim):
        """
        Des.
        """

        return self.dataContainer[iDim]

    def getDataArray(self, iDim):
        """
        Des.
        """

        if self.mpiComm != None:
            return self.dataContainer[iDim].getArray()
        else:
            return self.dataContainer[iDim]

    def assemble(self):
        """
        Des.
        """

        if self.mpiComm != None:
            for iDim in range(self.nDim):
                self.dataContainer[iDim].assemblyBegin()
                self.dataContainer[iDim].assemblyEnd()

    def norm(self):
        """
        Des.
        """

        normList = []

        for iDim in range(self.nDim):
            if self.mpiComm != None:
                val_norm = self.dataContainer[iDim].norm()
            else:
                val_norm = np.linalg.norm(self.dataContainer[iDim])
            normList.append(val_norm)

        return normList

    def sum(self):
        """
        Des.
        """

        sumList = []

        for iDim in range(self.nDim):
            val_sum = self.dataContainer[iDim].sum()
            sumList.append(val_sum)

        return sumList

    def copy(self):
        """
        Des.
        """

        newData = FlexInterfaceData(self.nPoint, self.nDim, self.mpiComm)
        for iDim in range(self.nDim):
            if self.mpiComm != None:
                self.dataContainer[iDim].copy(newData.dataContainer[iDim])
            else:
                np.copyto(newData.dataContainer[iDim], self.dataContainer[iDim])

        return newData

    def dot(self, data):
        """
        Des.
        """

        dotList = []

        if type(data) != type(self):
            raise TypeError("argument is not of type  FlexInterfaceData !")

        if self.nDim != data.nDim:
            raise IndexError("argument has not the same dimension !")
        else:
            for iDim in range(self.nDim):
                val_dot = self.dataContainer[iDim].dot(data.dataContainer[iDim])
                dotList.append(val_dot)

        return dotList

    def __setitem__(self, index, values):
        """
        Des.
        """

        if type(values) != list:
            raise TypeError("FlexInterfaceData.__setitem__ needs list as argument !")

        if len(values) != self.nDim:
            raise IndexError("Length of values does not match nDim !")
        else:
            for iDim in range(self.nDim):
                self.dataContainer[iDim][index] = values[iDim]

    #def __getitem__(self, index):
    #    """
    #    This is the implementation of the old InterfaceData class. Not reimplemented here because
    #    it is not really useful.
    #    """

    #    if self.mpiComm != None:
    #        send = None
    #        rcv = None
    #        for iProc in range(self.mpiSize):
    #            start, stop = self.indexing[iProc]
    #            if index in range(start, stop):
    #                sender = iProc
    #                break
    #        mpiBarrier(self.mpiComm)
    #        if self.myid == sender:
    #            send = (float(self.data_X[index]), float(self.data_Y[index]), float(self.data_Z[index]))
    #        rcv = self.mpiComm.bcast(send, sender)
    #        return rcv
    #    else:
    #        return (float(self.data_X[index]), float(self.data_Y[index]), float(self.data_Z[index]))

    def __add__(self, dataToAdd):
        """
        Des.
        """

        newData = FlexInterfaceData(self.nPoint, self.nDim, self.mpiComm)
        if type(dataToAdd) == type(self):
            if self.nDim != dataToAdd.nDim:
                raise IndexError("Dimensions do not match for + operator !")
            for iDim in range(self.nDim):
                newData.dataContainer[iDim] = self.dataContainer[iDim] + dataToAdd.dataContainer[iDim]
        elif type(dataToAdd) == float:
            for iDim in range(self.nDim):
                newData.dataContainer[iDim] = self.dataContainer[iDim] + dataToAdd
        elif type(dataToAdd) == int:
            for iDim in range(self.nDim):
                newData.dataContainer[iDim] = self.dataContainer[iDim] + float(dataToAdd)

        return newData

    def __radd__(self, dataToAdd):
        """
        Des.
        """

        newData = self + dataToAdd

        return newData

    def __iadd__(self, dataToAdd):
        """
        Des.
        """

        if type(dataToAdd) == type(self):
            if self.nDim != dataToAdd.nDim:
                raise IndexError("Dimensions do not match for += operator !")
            for iDim in range(self.nDim):
                self.dataContainer[iDim] += dataToAdd.dataContainer[iDim]
        elif type(dataToAdd) == float:
            for iDim in range(self.nDim):
                self.dataContainer[iDim] += dataToAdd
        elif type(dataToAdd) == int:
            for iDim in range(self.nDim):
                self.dataContainer[iDim] += float(dataToAdd)

        return self

    def __sub__(self, dataToSub):
        """
        Des.
        """

        newData = FlexInterfaceData(self.nPoint, self.nDim, self.mpiComm)
        if type(dataToSub) == type(self):
            if self.nDim != dataToSub.nDim:
                raise IndexError("Dimensions do not match for - operator !")
            for iDim in range(self.nDim):
                newData.dataContainer[iDim] = self.dataContainer[iDim] - dataToSub.dataContainer[iDim]
        elif type(dataToSub) == float:
            for iDim in range(self.nDim):
                newData.dataContainer[iDim] = self.dataContainer[iDim] - dataToSub
        elif type(dataToSub) == int:
            for iDim in range(self.nDim):
                newData.dataContainer[iDim] = self.dataContainer[iDim] - float(dataToSub)

        return newData

    def __rsub__(self, data):
        """
        Des.
        """

        newData = -1*self + data

        return newData

    def __isub__(self,dataToSub):
        """
        Des.
        """

        if type(dataToSub) == type(self):
            if self.nDim != dataToAdd.nDim:
                raise IndexError("Dimensions do not match for -= operator !")
            for iDim in range(self.nDim):
                self.dataContainer[iDim] -= dataToSub.dataContainer[iDim]
        elif type(dataToSub) == float:
            for iDim in range(self.nDim):
                self.dataContainer[iDim] -= dataToSub
        elif type(dataToSub) == int:
            for iDim in range(self.nDim):
                self.dataContainer[iDim] -= float(dataToSub)

        return self

    def __mul__(self, mulVal):
        """
        Des
        """

        newData = FlexInterfaceData(self.nPoint, self.nDim, self.mpiComm)
        for iDim in range(newData.nDim):
            newData.dataContainer[iDim] = self.dataContainer[iDim]*mulVal

        return newData

    def __rmul__(self, mulVal):
        """
        Des.
        """

        newData = self*mulVal

        return newData

    def __imul__(self, mulVal):
        """
        Des.
        """

        for iDim in range(self.nDim):
            self.dataContainer[iDim] *= mulVal

        return self
           

# ----------------------------------------------------------------------
#    InterfaceMatrix class
# ----------------------------------------------------------------------

class InterfaceMatrix:
    """
    Des.
    """

    def __init__(self, sizes, mpiComm=None):
        """
        Des.
        """

        self.sizes = sizes
        self.mpiComm = mpiComm

        if self.mpiComm != None:
            from petsc4py import PETSc
            self.H = PETSc.Mat().create(self.mpiComm)
            self.H.setType('mpiaij')
            self.H.setSizes(sizes)
            self.H.setUp()
        else:
            self.H = np.zeros(sizes, dtype=float)

    def setValue(self,(iGlobalIndex, jGlobalIndex), value):
        """
        Des.
        """

        if self.mpiComm != None:
            self.H.setValue(iGlobalIndex,jGlobalIndex,value)
        else:
            self.H[iGlobalIndex][jGlobalIndex] = value

    def assemble(self):
        """
        Des.
        """

        if self.mpiComm != None:
            self.H.assemblyBegin()
            self.H.assemblyEnd()

    def mult(self, Vec , VecOut):
        """
        Des.
        """

        if self.mpiComm != None:
            self.H.mult(Vec, VecOut)
        else:
            np.dot(self.H, Vec, VecOut)

    def getMat(self):
        """
        Des.
        """

        return self.H
                 
#class InterfaceDisplacement(InterfaceData)
#class InterfaceLoad(InterfaceData)

# ----------------------------------------------------------------------
#  Linear solver class
# ----------------------------------------------------------------------

class LinearSolver:
    """
    Description.
    Designed to be used with InterfaceData and InterfaceMatrix classes.
    """

    def __init__(self, MatrixOperator, mpiComm=None):
        """
        Description.
        MatrixOperator is of type InterfaceMatrix
        """
        self.mpiComm = mpiComm

        if mpiComm != None:
            from petsc4py import PETSc
            self.KSP_solver = PETSc.KSP().create(mpiComm)
            self.KSP_solver.setType('fgmres')
            self.KSP_solver.getPC().setType('jacobi')
            self.KSP_solver.setOperators(MatrixOperator.getMat())
            self.KSP_solver.setFromOptions()
            self.KSP_solver.setInitialGuessNonzero(True)
        else:
            self.LinOperator = MatrixOperator.getMat()
        

    def solve(self, VecX, VecB):
        """
        Des.
        """

        if self.mpiComm != None:
            self.KSP_solver.solve(VecX, VecB)
        else:
            VecB = splinalg.spsolve(self.LinOperator, VecX)

# ----------------------------------------------------------------------
#  Generic solid solver class
# ----------------------------------------------------------------------

class SolidSolver:
    """
    Des.
    """
    
    def __init__(self):
        
        self.haloNodeList = {}
        
        # --- Create the array for external communication (displacement, velocity and velocity at the previous time step) --- #
        self.nodalDisp_X = np.zeros(self.nPhysicalNodes)
        self.nodalDisp_Y = np.zeros(self.nPhysicalNodes)
        self.nodalDisp_Z = np.zeros(self.nPhysicalNodes)
        self.nodalVel_X = np.zeros(self.nPhysicalNodes)
        self.nodalVel_Y = np.zeros(self.nPhysicalNodes)
        self.nodalVel_Z = np.zeros(self.nPhysicalNodes)
        self.nodalVel_XNm1 = np.zeros(self.nPhysicalNodes)
        self.nodalVel_YNm1 = np.zeros(self.nPhysicalNodes)
        self.nodalVel_ZNm1 = np.zeros(self.nPhysicalNodes)

        # --- Same for thermal coupling (heat fluxes and temperature) ---
        self.nodalHeatFlux_X = np.zeros(self.nPhysicalNodes)
        self.nodalHeatFlux_Y = np.zeros(self.nPhysicalNodes)
        self.nodalHeatFlux_Z = np.zeros(self.nPhysicalNodes)
        self.nodalTemperature = np.zeros(self.nPhysicalNodes)
    
    def setInitialDisplacements(self):
        return
    
    def preprocessTimeIter(self, timeIter):
        return
    
    def run(self):
        return
    
    def __setCurrentState(self):
        return
    
    def getNodalDisplacements(self):
        """
        Des.
        """
        
        return (self.nodalDisp_X, self.nodalDisp_Y, self.nodalDisp_Z)

    def getNodalHeatFluxes(self):
        """
        Des.
        """

        return (self.nodalHeatFlux_X, self.nodalHeatFlux_Y, self.nodalHeatFlux_Z)

    def getNodalTemperatures(self):
        """
        Des.
        """

        return self.nodalTemperature
    
    def getNodalInitialPositions(self):
        return

    def getNodalVelocity(self):
        """
        des.
        """
        
        return (self.nodalVel_X, self.nodalVel_Y, self.nodalVel_Z)
    
    def getNodalVelocityNm1(self):
        """
        Des.
        """
        
        return (self.nodalVel_XNm1, self.nodalVel_YNm1, self.nodalVel_ZNm1)
    
    def getNodalIndex(self, iVertex):
        return
    
    def fakeFluidSolver(self, time):
        return
    
    def applyNodalLoads(self, load_X, load_Y, load_Z, time):
        return

    def applyNodalTemperatures(self, Temperature, time):
        return

    def applyNodalNormalHeatFluxes(self, NormalHeatFlux, val_time):
        return

    def applyNodalHeatFluxes(self, HeatFlux_X, HeatFlux_Y, HeatFlux_Z, time):
        return
    
    def update(self):
        
        self.nodalVel_XNm1 = self.nodalVel_X.copy()
        self.nodalVel_YNm1 = self.nodalVel_Y.copy()
        self.nodalVel_ZNm1 = self.nodalVel_Z.copy()

    def bgsUpdate(self):
        return
    
    def save(self):
        return
    
    def initRealTimeData(self):
        return
    
    def saveRealTimeData(self, time, nFSIIter):
        return

    def printRealTimeData(self, time, nFSIIter):
        return
    
    def remeshing(self):
        return
    
    def exit(self):
        return

# ----------------------------------------------------------------------
#  Generic fluid solver class
# ----------------------------------------------------------------------

class FluidSolver:
    """
    Des.
    """
    
    def __init__(self):
        
        self.haloNodeList = {}
        
        self.nodalLoad_X = np.zeros((self.nPhysicalNodes))
        self.nodalLoad_Y = np.zeros((self.nPhysicalNodes))
        self.nodalLoad_Z = np.zeros((self.nPhysicalNodes))

        self.nodalTemperature = np.zeros((self.nPhysicalNodes))

        self.nodalNormalHeatFlux = np.zeros(self.nPhysicalNodes)

        self.nodalHeatFlux_X = np.zeros((self.nPhysicalNodes))
        self.nodalHeatFlux_Y = np.zeros((self.nPhysicalNodes))
        self.nodalHeatFlux_Z = np.zeros((self.nPhysicalNodes))

        self.QWallInit = 0
        self.TWallInit = 288.0
    
    def setInitialMeshDeformation(self):
        return

    def setInitialInterfaceHeatFlux(self):
        return

    def setInitialInterfaceTemperature(self):
        return
    
    def preprocessTimeIter(self, timeIter):
        return
    
    def run(self):
        return
    
    def getNodalIndex(self, iVertex):
        return
    
    def fakeSolidSolver(self, time):
        return
    
    def getNodalLoads(self):
        
        return (self.nodalLoad_X, self.nodalLoad_Y, self.nodalLoad_Z)

    def getNodalInitialPositions(self):
        return

    def getNodalTemperatures(self):
        return self.nodalTemperature

    def getNodalNormalHeatFlux(self):
        return self.nodalNormalHeatFlux

    def getNodalHeatFluxes(self):
        return (self.nodalHeatFlux_X, self.nodalHeatFlux_Y, self.nodalHeatFlux_Z)
    
    def applyNodalDisplacements(self, dx, dy, dz, dx_nM1, dy_nM1, dz_nM1, haloNodesDisplacements,time):
        return

    def applyNodalHeatFluxes(self, HF_X, HF_Y, HF_Z, time):
        return

    def applyNodalTemperatures(self, Temperature, time):
        return
    
    def update(self, dt):
        return

    def bgsUpdate(self):
        return
    
    def save(self, nt):
        return
    
    def initRealTimeData(self):
        return
    
    def saveRealTimeData(self, time, nFSIIter):
        return

    def printRealTimeData(self, time, nFSIIter):
        return
    
    def remeshing(self):
        return
    
    def meshUpdate(self, nt):
        return
    
    def exit(self):
        return

# ----------------------------------------------------------------------
#    Criterion class
# ----------------------------------------------------------------------

class Criterion:
    """
    Description
    """
    
    def __init__(self, tolerance, heatFluxTolerance = 1e12):
        """
        Description.
        """
        
        self.tol = tolerance
        self.epsilon = 0

        self.tolHeatFlux = heatFluxTolerance
        self.epsilonHeatFlux = 0
    
    def isVerified(self, epsilon, epsilonHeatFlux=0):
        """
        Description.
        """

        verifList = [False, False]

        if epsilon < self.tol:
            verifList[0] = True

        if epsilonHeatFlux < self.tolHeatFlux:
            verifList[1] = True

        if False in verifList:
            return False
        else:
            return True

class DispNormCriterion(Criterion):
    """
    Description.
    """

    def __init__(self, tolerance, heatFluxTolerance = 1e12):
        """
        Description.
        """

        Criterion.__init__(self, tolerance, heatFluxTolerance)

    def update(self, res):
        """
        Des.
        """
        
        normX, normY, normZ = res.norm()
        
        norm = sqrt(normX**2 + normY**2 + normZ**2)
        
        self.epsilon = norm
        
        return self.epsilon

    def updateHeatFlux(self, resHeatFlux):
        """
        Des.
        """

        if resHeatFlux != None:
            #normX, normY, normZ = resHeatFlux.norm()
            normList = resHeatFlux.norm()
            normSquare = 0.0
            for ii in range(len(normList)):
                normSquare += normList[ii]**2

            #norm = sqrt(normX**2 + normY**2 + normZ**2)
            norm = sqrt(normSquare)
        else:
            norm = 1.0

        self.epsilonHeatFlux = norm

        return self.epsilonHeatFlux

# ----------------------------------------------------------------------
#    Manager class
# ----------------------------------------------------------------------

class Manager:
    """
    Description.
    """

    def __init__(self, FluidSolver, SolidSolver, nDim, computationType='steady', mpiComm=None):
        """
        Description.
        """
        
        mpiPrint('\n***************************** Initializing FSI interface *****************************', mpiComm)
        
        if mpiComm != None:
            self.mpiComm = mpiComm
            myid = mpiComm.Get_rank()
            mpiSize = mpiComm.Get_size()
        else:
            self.mpiComm = None
            myid = 0
            mpiSize = 1

        # --- Initialize all the parameters --- #
        self.nDim = nDim
        self.computationType = computationType
        self.withCht = False

        self.haveFluidSolver = False
        self.nLocalFluidInterfaceNodes = 0
        self.nLocalFluidInterfacePhysicalNodes = 0
        self.haveFluidInterface = False
        self.fluidHaloNodesList = {}
        self.fluidGlobalIndexRange = {}
        self.fluidIndexing = {}


        self.haveSolidSolver = False
        self.nLocalSolidInterfaceNodes = 0
        self.nLocalSolidInterfacePhysicalNodes = 0
        self.haveSolidInterface = False
        self.solidHaloNodesList = {}
        self.solidGlobalIndexRange = {}
        self.solidIndexing = {}

        # --- Identify the fluid and solid interfaces and store the number of nodes on both sides (and for each partition) ---

        if FluidSolver != None:
            print('Fluid solver is initialized on process {}'.format(myid))
            self.haveFluidSolver = True
            self.nLocalFluidInterfaceNodes = FluidSolver.nNodes
            if self.nLocalFluidInterfaceNodes != 0:
                self.haveFluidInterface = True
                print('Number of interface fluid nodes (halo nodes included) on proccess {} : {}'.format(myid,self.nLocalFluidInterfaceNodes))
        else:
            pass

        if SolidSolver != None:
            print('Solid solver is initialized on process {}'.format(myid))
            self.haveSolidSolver = True
            self.nLocalSolidInterfaceNodes = SolidSolver.nNodes
            if self.nLocalSolidInterfaceNodes != 0:
                self.haveSolidInterface = True
                print('Number of interface solid nodes (halo nodes included) on proccess {} : {}'.format(myid,self.nLocalSolidInterfaceNodes))
        else:
            pass

        # --- Exchange information about processors on which the solvers are defined and where the interface nodes are lying --- #
        if self.mpiComm != None:
            if self.haveFluidSolver == True:
                sendBufFluid = myid
            else:
                sendBufFluid = -1
            if self.haveSolidSolver == True:
                sendBufSolid = myid
            else:
                sendBufSolid = -1
            if self.haveFluidInterface == True:
                sendBufFluidInterface = myid
            else:
                sendBufFluidInterface = -1
            if self.haveSolidInterface == True:
                sendBufSolidInterface = myid
            else :
                sendBufSolidInterface = -1
            rcvBufFluid = mpiAllGather(mpiComm, sendBufFluid)
            rcvBufSolid = mpiAllGather(mpiComm, sendBufSolid)
            rcvBufFluidInterface = mpiAllGather(mpiComm, sendBufFluidInterface)
            rcvBufSolidInterface = mpiAllGather(mpiComm, sendBufSolidInterface)
            self.fluidSolverProcessors = rcvBufFluid[rcvBufFluid != -1]
            self.solidSolverProcessors = rcvBufSolid[rcvBufSolid != -1]
            self.fluidInterfaceProcessors = rcvBufFluidInterface[rcvBufFluidInterface != -1]
            self.solidInterfaceProcessors = rcvBufSolidInterface[rcvBufSolidInterface != -1]
        else:
            self.fluidSolverProcessors = np.zeros(1, dtype=int)
            self.solidSolverProcessors = np.zeros(1, dtype=int)
            self.fluidInterfaceProcessors = np.zeros(1, dtype=int)
            self.solidInterfaceProcessors = np.zeros(1, dtype=int)
        mpiBarrier(mpiComm)

        # --- Get the list of the halo nodes on the f/s interface --- #
        self.fluidHaloNodesList = FluidSolver.haloNodeList
        if myid in self.solidSolverProcessors:
            self.solidHaloNodesList = SolidSolver.haloNodeList
        if self.mpiComm != None:
            self.fluidHaloNodesList = self.mpiComm.allgather(self.fluidHaloNodesList)
            self.solidHaloNodesList = self.mpiComm.allgather(self.solidHaloNodesList)
        else:
            self.fluidHaloNodesList = [{}]
            self.solidHaloNodesList = [{}]

        # --- Get the number of physical (= not halo) nodes on the f/s interface --- #
        self.nLocalFluidInterfacePhysicalNodes = FluidSolver.nPhysicalNodes
        if myid in self.solidSolverProcessors:
            self.nLocalSolidInterfacePhysicalNodes = SolidSolver.nPhysicalNodes

        # --- Calculate the total (sum over all partitions) number of nodes at the f/s interface --- #
        self.nFluidInterfaceNodes = mpiAllReduce(mpiComm, self.nLocalFluidInterfaceNodes)
        self.nFluidInterfacePhysicalNodes = mpiAllReduce(mpiComm, self.nLocalFluidInterfacePhysicalNodes)
        self.nSolidInterfaceNodes = mpiAllReduce(mpiComm, self.nLocalSolidInterfaceNodes)
        self.nSolidInterfacePhysicalNodes = mpiAllReduce(mpiComm, self.nLocalSolidInterfacePhysicalNodes)
        mpiPrint('Total number of fluid interface nodes (halo nodes included) : {}'.format(self.nFluidInterfaceNodes), mpiComm)
        mpiPrint('Total number of solid interface nodes (halo nodes included) : {}'.format(self.nSolidInterfaceNodes), mpiComm)
        mpiPrint('Total number of fluid interface nodes : {}'.format(self.nFluidInterfacePhysicalNodes), mpiComm)
        mpiPrint('Total number of solid interface nodes : {}'.format(self.nSolidInterfacePhysicalNodes), mpiComm)

        # --- Store the number of physical interface nodes on each processor and allgather the information --- #
        self.fluidPhysicalInterfaceNodesDistribution = np.zeros(mpiSize, dtype=int)
        self.solidPhysicalInterfaceNodesDistribution = np.zeros(mpiSize, dtype=int)
        if self.mpiComm != None:
            self.fluidPhysicalInterfaceNodesDistribution = mpiAllGather(self.mpiComm, self.nLocalFluidInterfacePhysicalNodes)
            self.solidPhysicalInterfaceNodesDistribution = mpiAllGather(self.mpiComm, self.nLocalSolidInterfacePhysicalNodes)
        else:
            self.fluidPhysicalInterfaceNodesDistribution[0] = self.nFluidInterfacePhysicalNodes
            self.solidPhysicalInterfaceNodesDistribution[0] = self.nSolidInterfacePhysicalNodes

        # --- Calculate and store the global indexing of interface physical nodes
        if self.mpiComm != None:
            self.fluidGlobalIndexRange = {}
            self.solidGlobalIndexRange = {}
            if myid in self.fluidInterfaceProcessors:
                globalIndexStart = 0
                for iProc in range(myid):
                    globalIndexStart += self.fluidPhysicalInterfaceNodesDistribution[iProc]
                globalIndexStop = globalIndexStart + self.nLocalFluidInterfacePhysicalNodes-1
            else:
                globalIndexStart = 0
                globalIndexStop = 0
            self.fluidGlobalIndexRange[myid] = [globalIndexStart,globalIndexStop]
            self.fluidGlobalIndexRange = self.mpiComm.allgather(self.fluidGlobalIndexRange)
            if myid in self.solidInterfaceProcessors:
                globalIndexStart = 0
                for jProc in range(myid):
                    globalIndexStart += self.solidPhysicalInterfaceNodesDistribution[jProc]
                globalIndexStop = globalIndexStart + self.nLocalSolidInterfaceNodes-1
            else:
                globalIndexStart = 0
                globalIndexStop = 0
            self.solidGlobalIndexRange[myid] = [globalIndexStart,globalIndexStop]
            self.solidGlobalIndexRange = self.mpiComm.allgather(self.solidGlobalIndexRange)
        else:
            temp = {}
            temp[0] = [0,self.nLocalFluidInterfacePhysicalNodes-1]
            self.fluidGlobalIndexRange = list()
            self.fluidGlobalIndexRange.append(temp)
            temp[0] = [0,self.nSolidInterfacePhysicalNodes-1]
            self.solidGlobalIndexRange = list()
            self.solidGlobalIndexRange.append(temp)

        # --- Map the FSI indexing with the solvers indexing --- #
        fluidIndexing_temp = {}
        localIndex = 0
        for iVertex in range(self.nLocalFluidInterfaceNodes):
            nodeIndex = FluidSolver.getNodalIndex(iVertex)
            if nodeIndex in self.fluidHaloNodesList[myid].keys():
                pass
            else:
                fluidIndexing_temp[nodeIndex] = self.getGlobalIndex('fluid', myid, localIndex)
                localIndex += 1

        solidIndexing_temp = {}
        localIndex = 0
        for jVertex in range(self.nLocalSolidInterfaceNodes):
            nodeIndex = SolidSolver.getNodalIndex(jVertex)
            if nodeIndex in self.solidHaloNodesList[myid].keys():
                pass
            else:
                solidIndexing_temp[nodeIndex] = self.getGlobalIndex('solid', myid, localIndex)
                localIndex += 1

        if self.mpiComm != None:
            fluidIndexing_temp = self.mpiComm.allgather(fluidIndexing_temp)
            solidIndexing_temp = self.mpiComm.allgather(solidIndexing_temp)
            for ii in range(len(solidIndexing_temp)):
                for key, value in solidIndexing_temp[ii].items():
                    self.solidIndexing[key] = value
            for ii in range(len(fluidIndexing_temp)):
                for key, value in fluidIndexing_temp[ii].items():
                    self.fluidIndexing[key] = value
        else:
            self.fluidIndexing = fluidIndexing_temp.copy()
            self.solidIndexing = solidIndexing_temp.copy()
        del fluidIndexing_temp, solidIndexing_temp
        
    def getGlobalIndex(self, domain, iProc, iLocalVertex):
        """
        Description.
        """

        if domain == 'fluid':
            globalStartIndex = self.fluidGlobalIndexRange[iProc][iProc][0]
        elif domain == 'solid':
            globalStartIndex = self.solidGlobalIndexRange[iProc][iProc][0]
        globalIndex = globalStartIndex + iLocalVertex

        return globalIndex

    def getNumberOfFluidInterfaceNodes(self):
        """
        Description.
        """

        return self.nFluidInterfacePhysicalNodes

    def getNumberOfLocalFluidInterfaceNodes(self):
        """
        Description.
        """

        return self.nLocalFluidInterfacePhysicalNodes

    def getNumberOfSolidInterfaceNodes(self):
        """
        Description.
        """

        return self.nSolidInterfacePhysicalNodes

    def getNumberOfLocalSolidInterfaceNodes(self):
        """
        Description.
        """

        return self.nLocalSolidInterfacePhysicalNodes

    def getSolidSolverProcessors(self):
        """
        Des.
        """

        return self.solidSolverProcessors

    def getSolidInterfaceProcessors(self):
        """
        Description.
        """

        return self.solidInterfaceProcessors

    def getFluidInterfaceProcessors(self):
        """
        Description.
        """

        return self.fluidInterfaceProcessors

    def getSolidPhysicalInterfaceNodesDistribution(self):
        """
        Des.
        """

        return self.solidPhysicalInterfaceNodesDistribution

    def getFluidPhysicalInterfaceNodesDistribution(self):
        """
        Des.
        """

        return self.fluidPhysicalInterfaceNodesDistribution

    def getSolidGlobalIndexRange(self):
        """
        Des.
        """

        return self.solidGlobalIndexRange

    def getFluidGlobalIndexRange(self):
        """
        Des.
        """

        return self.fluidGlobalIndexRange

    def getFluidHaloNodesList(self):
        """
        Des.
        """

        return self.fluidHaloNodesList

    def getSolidHaloNodesList(self):
        """
        Des.
        """

        return self.solidHaloNodesList

    def getFluidIndexing(self):
        """
        Des.
        """

        return self.fluidIndexing

    def getSolidIndexing(self):
        """
        Des.
        """

        return self.solidIndexing

    def getnDim(self):
        """
        Des.
        """
 
        return self.nDim

    def getComputationType(self):
        """
        des.
        """

        return self.computationType

    def getMPIComm(self):
        """
        Description.
        """

        return self.mpiComm

# ----------------------------------------------------------------------
#    Interpolator class
# ----------------------------------------------------------------------

class InterfaceInterpolator:
    """
    Description.
    """

    def __init__(self, Manager, FluidSolver, SolidSolver, mpiComm = None):
        """
        Description.
        """
        
        mpiPrint('\n***************************** Initializing FSI interpolator *****************************', mpiComm)
        
        self.manager = Manager
        self.SolidSolver = SolidSolver
        self.FluidSolver = FluidSolver

        self.nf = 0
        self.ns = 0
        self.nf_loc = 0
        self.ns_loc = 0
        self.nDim = 0

        self.d = 0

        self.chtTransferMethod = None
        #FFTB = Flux Forward Temperature Back
        #TFFB = Temperature Forward Flux Back
        #hFTB = Heat Transfer Coefficient Forward Temperature Back
        #hFFB = Heat Transfer Coefficient Forward Flux Back
        self.heatTransferCoeff = None

        self.mpiComm = mpiComm
        
        if self.mpiComm != None:
            self.myid = self.mpiComm.Get_rank()
            self.mpiSize = self.mpiComm.Get_size()
        else:
            self.myid = 0
            self.mpiSize = 1
 
    def generateInterfaceData(self):
        """
        Description.
        """

        self.solidInterfaceDisplacement = FlexInterfaceData(self.ns + self.d, 3, self.mpiComm)
        self.fluidInterfaceDisplacement = FlexInterfaceData(self.nf, 3, self.mpiComm)
        self.solidInterfaceLoads = FlexInterfaceData(self.ns + self.d, 3, self.mpiComm)
        self.fluidInterfaceLoads = FlexInterfaceData(self.nf, 3, self.mpiComm)

        self.solidInterfaceHeatFlux = None
        self.fluidInterfaceHeatFlux = None
        self.solidInterfaceTemperature = None
        self.fluidInterfaceTemperature = None
        self.fluidInterfaceNormalHeatFlux = None
        self.solidInterfaceNormalHeatFlux = None
        self.fluidInterfaceRobinTemperature = None
        self.solidInterfaceRobinTemperature = None
        if self.manager.withCht :
            self.solidInterfaceHeatFlux = FlexInterfaceData(self.ns + self.d, 3, self.mpiComm)
            self.fluidInterfaceHeatFlux = FlexInterfaceData(self.nf, 3, self.mpiComm)
            self.solidInterfaceTemperature = FlexInterfaceData(self.ns + self.d, 1, self.mpiComm)
            self.fluidInterfaceTemperature = FlexInterfaceData(self.nf, 1, self.mpiComm)
            self.fluidInterfaceNormalHeatFlux = FlexInterfaceData(self.nf, 1, self.mpiComm)
            self.solidInterfaceNormalHeatFlux = FlexInterfaceData(self.ns, 1, self.mpiComm)
            #if self.chtTransferMethod == 'hFFB' or self.chtTransferMethod == 'hFTB':
            self.fluidInterfaceRobinTemperature = FlexInterfaceData(self.nf, 1, self.mpiComm)
            self.solidInterfaceRobinTemperature = FlexInterfaceData(self.ns, 1, self.mpiComm)


    def createKDTree(self, posX_array, posY_array, posZ_array):
        """
        Description.
        """

        KDTree = spatial.KDTree(zip(posX_array, posY_array, posZ_array))

        return KDTree

    def distance(self, NodeA, NodeB):
        """
        Des.
        """

        return sqrt((NodeB[0] - NodeA[0])**2 + (NodeB[1] - NodeA[1])**2 + (NodeB[2] - NodeA[2])**2)

    def checkConservation(self):
        """
        Des.
        """

        WSX, WSY, WSZ = self.solidInterfaceLoads.dot(self.solidInterfaceDisplacement)

        WFX, WFY, WFZ = self.fluidInterfaceLoads.dot(self.fluidInterfaceDisplacement)

        mpiPrint("Checking f/s interface conservation...", self.mpiComm)
        mpiPrint('Solid side (Wx, Wy, Wz) = ({}, {}, {})'.format(WSX, WSY, WSZ), self.mpiComm)                
        mpiPrint('Fluid side (Wx, Wy, Wz) = ({}, {}, {})'.format(WFX, WFY, WFZ), self.mpiComm)

    def checkTotalLoad(self):
        """
        Des.
        """

        FX, FY, FZ = self.solidInterfaceLoads.sum()

        FFX, FFY, FFZ = self.fluidInterfaceLoads.sum()

        mpiPrint("Checking f/s interface total force...", self.mpiComm)
        mpiPrint('Solid side (Fx, Fy, Fz) = ({}, {}, {})'.format(FX, FY, FZ), self.mpiComm)                
        mpiPrint('Fluid side (Fx, Fy, Fz) = ({}, {}, {})'.format(FFX, FFY, FFZ), self.mpiComm)

    def getDisplacementFromSolidSolver(self):
        """
        Des.
        """
        
        if self.myid in self.manager.getSolidInterfaceProcessors():
            localSolidInterfaceDisp_X, localSolidInterfaceDisp_Y, localSolidInterfaceDisp_Z = self.SolidSolver.getNodalDisplacements()
            for iVertex in range(self.ns_loc):
                iGlobalVertex = self.manager.getGlobalIndex('solid', self.myid, iVertex)
                self.solidInterfaceDisplacement[iGlobalVertex] = [localSolidInterfaceDisp_X[iVertex], localSolidInterfaceDisp_Y[iVertex], localSolidInterfaceDisp_Z[iVertex]]

        self.solidInterfaceDisplacement.assemble()

    def getHeatFluxFromSolidSolver(self):
        """
        Des.
        """

        if self.myid in self.manager.getSolidInterfaceProcessors():
            localSolidInterfaceHeatFlux_X, localSolidInterfaceHeatFlux_Y, localSolidInterfaceHeatFlux_Z = self.SolidSolver.getNodalHeatFluxes()
            for iVertex in range(self.ns_loc):
                iGlobalVertex = self.manager.getGlobalIndex('solid', self.myid, iVertex)
                self.solidInterfaceHeatFlux[iGlobalVertex] = [localSolidInterfaceHeatFlux_X[iVertex], localSolidInterfaceHeatFlux_Y[iVertex], localSolidInterfaceHeatFlux_Z[iVertex]]

        self.solidInterfaceHeatFlux.assemble()

    def getLoadsFromFluidSolver(self):
        """
        Des.
        """

        if self.myid in self.manager.getFluidInterfaceProcessors():
            localFluidInterfaceLoad_X, localFluidInterfaceLoad_Y, localFluidInterfaceLoad_Z = self.FluidSolver.getNodalLoads()
            for iVertex in range(self.nf_loc):
                iGlobalVertex = self.manager.getGlobalIndex('fluid', self.myid, iVertex)
                self.fluidInterfaceLoads[iGlobalVertex] = [localFluidInterfaceLoad_X[iVertex], localFluidInterfaceLoad_Y[iVertex], localFluidInterfaceLoad_Z[iVertex]]

        self.fluidInterfaceLoads.assemble()

    def getTemperatureFromFluidSolver(self):
        """
        Des.
        """

        if self.myid in self.manager.getFluidInterfaceProcessors():
            localFluidInterfaceTemperature = self.FluidSolver.getNodalTemperatures()
            for iVertex in range(self.nf_loc):
                iGlobalVertex = self.manager.getGlobalIndex('fluid', self.myid, iVertex)
                self.fluidInterfaceTemperature[iGlobalVertex] = [localFluidInterfaceTemperature[iVertex]]

        self.fluidInterfaceTemperature.assemble()

    def getRobinTemperatureFromFluidSolver(self):
        """
        Des.
        """

        if self.myid in self.manager.getFluidInterfaceProcessors():
            localFluidInterfaceNormalHeatFlux = self.FluidSolver.getNodalNormalHeatFlux()
            localFluidInterfaceTemperature = self.FluidSolver.getNodalTemperatures()
            localFluidInterfaceRobinTemperature = localFluidInterfaceTemperature - (localFluidInterfaceNormalHeatFlux/self.heatTransferCoeff)
            for iVertex in range(self.nf_loc):
                iGlobalVertex = self.manager.getGlobalIndex('fluid', self.myid, iVertex)
                self.fluidInterfaceRobinTemperature[iGlobalVertex] = [localFluidInterfaceRobinTemperature[iVertex]]

        self.fluidInterfaceRobinTemperature.assemble()

    def getHeatFluxFromFluidSolver(self):
        """
        Des.
        """

        if self.myid in self.manager.getFluidInterfaceProcessors():
            localFluidInterfaceHeatFlux_X, localFluidInterfaceHeatFlux_Y, localFluidInterfaceHeatFlux_Z = self.FluidSolver.getNodalHeatFluxes()
            localFluidInterfaceNormalHeatFlux = self.FluidSolver.getNodalNormalHeatFlux()
            print localFluidInterfaceNormalHeatFlux
            for iVertex in range(self.nf_loc):
                iGlobalVertex = self.manager.getGlobalIndex('fluid', self.myid, iVertex)
                self.fluidInterfaceHeatFlux[iGlobalVertex] = [localFluidInterfaceHeatFlux_X[iVertex], localFluidInterfaceHeatFlux_Y[iVertex], localFluidInterfaceHeatFlux_Z[iVertex]]
                self.fluidInterfaceNormalHeatFlux[iGlobalVertex] = [localFluidInterfaceNormalHeatFlux[iVertex]]

        self.fluidInterfaceHeatFlux.assemble()
        self.fluidInterfaceNormalHeatFlux.assemble()

    def redistributeDataToFluidSolver(self, fluidInterfaceData):
        """
        Description
        """

        localFluidInterfaceData_array = None
        haloNodesData = {}

        if self.mpiComm != None:
            localSize = fluidInterfaceData.getDataArray(0).shape[0]
            fluidInterfaceData_array_recon = []
            for iDim in range(fluidInterfaceData.nDim):
                array_recon = mpiGatherv(fluidInterfaceData.getDataArray(iDim), localSize, self.nf, self.mpiComm, 0)
                fluidInterfaceData_array_recon.append(array_recon)
            haloNodesData = {}
            haloNodesData_bis = {}
            if self.myid == 0:
                for iProc in self.manager.getFluidInterfaceProcessors():
                    fluidPhysicalInterfaceNodesDistribution = self.manager.getFluidPhysicalInterfaceNodesDistribution()
                    fluidGlobalIndexRange = self.manager.getFluidGlobalIndexRange()
                    sendBuff = []
                    for iDim in range(fluidInterfaceData.nDim):
                        sendBuff_i = np.zeros(fluidPhysicalInterfaceNodesDistribution[iProc])
                        sendBuff.append(sendBuff_i)
                    globalIndex = fluidGlobalIndexRange[iProc][iProc][0]
                    sendBuffHalo = {}
                    for iVertex in range(fluidPhysicalInterfaceNodesDistribution[iProc]):
                        for iDim in range(fluidInterfaceData.nDim):
                            sendBuff[iDim][iVertex] = fluidInterfaceData_array_recon[iDim][globalIndex]
                        globalIndex += 1
                    fluidHaloNodesList = self.manager.getFluidHaloNodesList()
                    fluidIndexing = self.manager.getFluidIndexing()
                    for key in fluidHaloNodesList[iProc].keys():
                        globalIndex = fluidIndexing[key]
                        sendBuffHalo[key] = []
                        for iDim in range(fluidInterfaceData.nDim):
                            sendBuffHalo[key].append(fluidInterfaceData_array_recon[iDim][globalIndex])
                    iTagSend = 1
                    for iDim in range(fluidInterfaceData.nDim):
                        self.mpiComm.Send(sendBuff[iDim], dest=iProc, tag = iTagSend)
                        iTagSend += 1
                    #self.mpiComm.send(sendBuffHalo, dest = iProc, tag=iTagSend)
                    sendBuffHalo_key = np.array(sendBuffHalo.keys())
                    sendBuffHalo_values = np.empty((sendBuffHalo_key.size, 3),dtype=float)
                    for ii in range(sendBuffHalo_key.size):
                        sendBuffHalo_values[ii] = np.array(sendBuffHalo[sendBuffHalo_key[ii]])
                    self.mpiComm.Send(np.array(sendBuffHalo_key.size), dest=iProc, tag=101)
                    self.mpiComm.Send(sendBuffHalo_key, dest=iProc, tag=102)
                    self.mpiComm.Send(sendBuffHalo_values, dest=iProc, tag=103)
            if self.myid in self.manager.getFluidInterfaceProcessors():
                localFluidInterfaceData_array = []
                iTagRec = 1
                for iDim in range(fluidInterfaceData.nDim):
                    local_array = np.zeros(self.nf_loc)
                    self.mpiComm.Recv(local_array, source=0, tag=iTagRec)
                    localFluidInterfaceData_array.append(local_array)
                    iTagRec += 1
                #haloNodesData = self.mpiComm.recv(source=0, tag=iTagRec)
                nHaloNodesRcv = np.empty(1, dtype=int)
                self.mpiComm.Recv(nHaloNodesRcv, source=0, tag=101)
                rcvBuffHalo_keyBuff = np.empty(nHaloNodesRcv[0], dtype=int)
                self.mpiComm.Recv(rcvBuffHalo_keyBuff, source=0, tag=102)
                rcvBuffHalo_values = np.empty((nHaloNodesRcv[0],3), dtype=float)
                self.mpiComm.Recv(rcvBuffHalo_values, source=0, tag=103)
                for ii in range(len(rcvBuffHalo_keyBuff)):
                    haloNodesData_bis[rcvBuffHalo_keyBuff[ii]] = list(rcvBuffHalo_values[ii])
                haloNodesData = haloNodesData_bis
                

        return (localFluidInterfaceData_array, haloNodesData)

    def redistributeDataToSolidSolver(self, solidInterfaceData):
        """
        Des.
        """

        localSolidInterfaceData_array = None
        haloNodesData = {}
        haloNodesData_bis = {}

        if self.mpiComm != None:
            localSize = solidInterfaceData.getDataArray(0).shape[0]
            solidInterfaceData_array_recon = []
            for iDim in range(solidInterfaceData.nDim):
                array_recon = mpiGatherv(solidInterfaceData.getDataArray(iDim), localSize, self.ns+self.d, self.mpiComm, 0)
                solidInterfaceData_array_recon.append(array_recon)
            haloNodesData = {}
            if self.myid == 0:
                for iProc in self.manager.getSolidInterfaceProcessors():
                    solidPhysicalInterfaceNodesDistribution = self.manager.getSolidPhysicalInterfaceNodesDistribution()
                    solidGlobalIndexRange = self.manager.getSolidGlobalIndexRange()
                    sendBuff = []
                    for iDim in range(solidInterfaceData.nDim):
                        sendBuff_i = np.zeros(solidPhysicalInterfaceNodesDistribution[iProc])
                        sendBuff.append(sendBuff_i)
                    globalIndex = solidGlobalIndexRange[iProc][iProc][0]
                    sendBuffHalo = {}
                    for iVertex in range(solidPhysicalInterfaceNodesDistribution[iProc]):
                        for iDim in range(solidInterfaceData.nDim):
                            sendBuff[iDim][iVertex] = solidInterfaceData_array_recon[iDim][globalIndex]
                        globalIndex += 1
                    solidHaloNodesList = self.manager.getSolidHaloNodesList()
                    solidIndexing = self.manager.getSolidIndexing()
                    for key in solidHaloNodesList[iProc].keys():
                        globalIndex = solidIndexing[key]
                        sendBuffHalo[key] = []
                        for iDim in range(solidInterfaceData.nDim):
                            sendBuffHalo[key].append(solidInterfaceData_array_recon[iDim][globalIndex])
                    iTagSend = 1
                    for iDim in range(solidInterfaceData.nDim):
                        self.mpiComm.Send(sendBuff[iDim], dest=iProc, tag = iTagSend)
                        iTagSend += 1
                    #self.mpiComm.send(sendBuffHalo, dest = iProc, tag=iTagSend)
                    sendBuffHalo_key = np.array(sendBuffHalo.keys())
                    sendBuffHalo_values = np.empty((sendBuffHalo_key.size, 3),dtype=float)
                    for ii in range(sendBuffHalo_key.size):
                        sendBuffHalo_values[ii] = np.array(sendBuffHalo[sendBuffHalo_key[ii]])
                    self.mpiComm.Send(np.array(sendBuffHalo_key.size), dest=iProc, tag=101)
                    self.mpiComm.Send(sendBuffHalo_key, dest=iProc, tag=102)
                    self.mpiComm.Send(sendBuffHalo_values, dest=iProc, tag=103)
            if self.myid in self.manager.getSolidInterfaceProcessors():
                localSolidInterfaceData_array = []
                iTagRec = 1
                for iDim in range(solidInterfaceData.nDim):
                    local_array = np.zeros(self.ns_loc)
                    self.mpiComm.Recv(local_array, source=0, tag = iTagRec)
                    localSolidInterfaceData_array.append(local_array)
                    iTagRec += 1
                #haloNodesData = self.mpiComm.recv(source=0, tag=iTagRec)
                nHaloNodesRcv = np.empty(1, dtype=int)
                self.mpiComm.Recv(nHaloNodesRcv, source=0, tag=101)
                rcvBuffHalo_keyBuff = np.empty(nHaloNodesRcv[0], dtype=int)
                self.mpiComm.Recv(rcvBuffHalo_keyBuff, source=0, tag=102)
                rcvBuffHalo_values = np.empty((nHaloNodesRcv[0],3), dtype=float)
                self.mpiComm.Recv(rcvBuffHalo_values, source=0, tag=103)
                for ii in range(len(rcvBuffHalo_keyBuff)):
                    haloNodesData_bis[rcvBuffHalo_keyBuff[ii]] = list(rcvBuffHalo_values[ii])
                haloNodesData = haloNodesData_bis

        return (localSolidInterfaceData_array, haloNodesData)

    def setLoadsToSolidSolver(self, time):
        """
        des.
        """

        FFX, FFY, FFZ = self.fluidInterfaceLoads.sum()

        
        FX = 0.
        FY = 0.
        FZ = 0.

        FXT = 0.
        FYT = 0.
        FZT = 0.

        if self.mpiComm != None:
            (localSolidLoads_array, haloNodesSolidLoads) = self.redistributeDataToSolidSolver(self.solidInterfaceLoads)
            if self.myid in self.manager.getSolidInterfaceProcessors():
                self.SolidSolver.applyNodalLoads(localSolidLoads_array[0], localSolidLoads_array[1], localSolidLoads_array[2], time)
                FX = localSolidLoads_array[0].sum()
                FY = localSolidLoads_array[1].sum()
                FZ = localSolidLoads_array[2].sum()
            FXT = mpiAllReduce(self.mpiComm, FX)
            FYT = mpiAllReduce(self.mpiComm, FY)
            FZT = mpiAllReduce(self.mpiComm, FZ)
        else:
            self.SolidSolver.applyNodalLoads(self.solidInterfaceLoads.getDataArray(0), self.solidInterfaceLoads.getDataArray(1), self.solidInterfaceLoads.getDataArray(2), time)
            FXT, FYT, FZT = self.solidInterfaceLoads.sum()

        mpiPrint("Checking f/s interface total force...", self.mpiComm)
        mpiPrint('Solid side (Fx, Fy, Fz) = ({}, {}, {})'.format(FXT, FYT, FZT), self.mpiComm)                
        mpiPrint('Fluid side (Fx, Fy, Fz) = ({}, {}, {})'.format(FFX, FFY, FFZ), self.mpiComm)

    def setDisplacementToFluidSolver(self, time):
        """
        Des.
        """

        self.checkConservation()

        if self.mpiComm != None:
            (localFluidInterfaceDisplacement, haloNodesDisplacements) = self.redistributeDataToFluidSolver(self.fluidInterfaceDisplacement)
            if self.myid in self.manager.getFluidInterfaceProcessors():
                self.FluidSolver.applyNodalDisplacements(localFluidInterfaceDisplacement[0], localFluidInterfaceDisplacement[1], localFluidInterfaceDisplacement[2], localFluidInterfaceDisplacement[0], localFluidInterfaceDisplacement[1], localFluidInterfaceDisplacement[2], haloNodesDisplacements, time)
        else:
            self.FluidSolver.applyNodalDisplacements(self.fluidInterfaceDisplacement.getDataArray(0), self.fluidInterfaceDisplacement.getDataArray(1), self.fluidInterfaceDisplacement.getDataArray(2), self.fluidInterfaceDisplacement.getDataArray(0), self.fluidInterfaceDisplacement.getDataArray(1), self.fluidInterfaceDisplacement.getDataArray(2), {}, time)

    def setHeatFluxToFluidSolver(self, time):
        """
        Description.
        """

        if self.mpiComm != None:
            (localFluidInterfaceHeatFlux, haloNodesHeatFlux) = self.redistributeDataToFluidSolver(self.fluidInterfaceHeatFlux)
            if self.myid in self.manager.getFluidInterfaceProcessors():
                self.FluidSolver.applyNodalHeatFluxes(localFluidInterfaceHeatFlux[0], localFluidInterfaceHeatFlux[1], localFluidInterfaceHeatFlux[2], time)
        else:
            self.FluidSolver.applyNodalHeatFluxes(self.fluidInterfaceHeatFlux.getDataArray(0), self.fluidInterfaceHeatFlux.getDataArray(1), self.fluidInterfaceHeatFlux.getDataArray(2), time)

    def setTemperatureToFluidSolver(self, time):
        """
        Des.
        """

        if self.mpiComm != None:
            (localFluidInterfaceTemperature, haloNodesTemperature) = self.redistributeDataToFluidSolver(self.fluidInterfaceTemperature)
            if self.myid in self.manager.getFluidInterfaceProcessors():
                self.FluidSolver.applyNodalTemperatures(localFluidInterfaceTemperature[0], time)
        else:
            self.FluidSolver.applyNodalTemperatures(self.fluidInterfaceTemperature.getDataArray(0), time)

    def setTemperatureToSolidSolver(self, time):
        """
        Description
        """

        if self.mpiComm != None:
            (localSolidInterfaceTemperature, haloNodesTemperature) = self.redistributeDataToSolidSolver(self.solidInterfaceTemperature)
            if self.myid in self.manager.getSolidInterfaceProcessors():
                self.SolidSolver.applyNodalTemperatures(localSolidInterfaceTemperature[0], time)
        else:
            self.SolidSolver.applyNodalTemperatures(self.solidInterfaceTemperature.getDataArray(0), time)

    def setRobinHeatFluxToSolidSolver(self, time):
        """
        Def
        """

        if self.mpiComm != None:
            (localSolidInterfaceRobinTemperature, haloNodesRobinTemperature) = self.redistributeDataToSolidSolver(self.solidInterfaceRobinTemperature)
            if self.myid in self.manager.getSolidInterfaceProcessors():
                localSolidInterfaceTemperature = self.SolidSolver.getNodalTemperatures()
                localSolidInterfaceRobinHeatFlux = self.heatTransferCoeff*(localSolidInterfaceTemperature-localSolidInterfaceRobinTemperature[0])
                self.SolidSolver.applyNodalNormalHeatFluxes(localSolidInterfaceRobinHeatFlux, time)
        else:
            localSolidInterfaceTemperature = self.SolidSolver.getNodalTemperatures()
            localSolidInterfaceRobinHeatFlux = self.heatTransferCoeff*(localSolidInterfaceTemperature-self.solidInterfaceRobinTemperature.getDataArray(0), time)
            self.SolidSolver.applyNodalNormalHeatFluxes(localSolidInterfaceRobinHeatFlux, time)

    def setHeatFluxToSolidSolver(self, time):
        """
        Des.
        """

        if self.mpiComm != None:
            (localSolidInterfaceNormalHeatFlux, haloNodesNormalHeatFlux) =  self.redistributeDataToSolidSolver(self.solidInterfaceNormalHeatFlux)
            if self.myid in self.manager.getSolidInterfaceProcessors():
                self.SolidSolver.applyNodalNormalHeatFluxes(localSolidInterfaceNormalHeatFlux[0], time)
        else:
            self.SolidSolver.applyNodalNormalHeatFluxes(self.solidInterfaceNormalHeatFlux.getDataArray(0), time)

    def getNs(self):
        """
        Des.
        """

        return self.ns

    def getNf(self):
        """
        Des.
        """

        return self.nf

    def getd(self):
        """
        Des.
        """

        return self.d

class MatchingMeshesInterpolator(InterfaceInterpolator):
    """
    Description.
    """

    def __init__(self, Manager, FluidSolver, SolidSolver, mpiComm = None):
        """
        Description
        """

        InterfaceInterpolator.__init__(self, Manager, FluidSolver, SolidSolver, mpiComm)

        mpiPrint('\nSetting matching meshes interpolator...', mpiComm)

        self.nf = self.manager.getNumberOfFluidInterfaceNodes()
        self.ns = self.manager.getNumberOfSolidInterfaceNodes()
        self.nf_loc = self.manager.getNumberOfLocalFluidInterfaceNodes()
        self.ns_loc = self.manager.getNumberOfLocalSolidInterfaceNodes()
        self.nDim = self.manager.getnDim()


        if self.nf != self.ns:
            raise Exception("Fluid and solid interface must have the same number of nodes for matching meshes ! ")

        self.generateInterfaceData()

        self.H = InterfaceMatrix((self.nf,self.ns), self.mpiComm)
        self.H_T = InterfaceMatrix((self.ns,self.nf), self.mpiComm)
                 
        solidInterfaceProcessors = self.manager.getSolidInterfaceProcessors()
        fluidInterfaceProcessors = self.manager.getFluidInterfaceProcessors()
        solidPhysicalInterfaceNodesDistribution = self.manager.getSolidPhysicalInterfaceNodesDistribution()

        mpiPrint('\nBuilding interpolation matrix...', mpiComm)

        if self.mpiComm != None:
            for iProc in solidInterfaceProcessors:
                if self.myid == iProc:
                    localSolidInterface_array_X, localSolidInterface_array_Y, localSolidInterface_array_Z = self.SolidSolver.getNodalInitialPositions()
                    for jProc in fluidInterfaceProcessors:
                        self.mpiComm.Send(localSolidInterface_array_X, dest=jProc, tag=1)
                        self.mpiComm.Send(localSolidInterface_array_Y, dest=jProc, tag=2)
                        self.mpiComm.Send(localSolidInterface_array_Z, dest=jProc, tag=3)
                if self.myid in fluidInterfaceProcessors:
                    sizeOfBuff = solidPhysicalInterfaceNodesDistribution[iProc]
                    solidInterfaceBuffRcv_X = np.zeros(sizeOfBuff)
                    solidInterfaceBuffRcv_Y = np.zeros(sizeOfBuff)
                    solidInterfaceBuffRcv_Z = np.zeros(sizeOfBuff)
                    self.mpiComm.Recv(solidInterfaceBuffRcv_X, iProc, tag=1)
                    self.mpiComm.Recv(solidInterfaceBuffRcv_Y, iProc, tag=2)
                    self.mpiComm.Recv(solidInterfaceBuffRcv_Z, iProc, tag=3)
                    self.fillMatrix(solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc)
        else:
            localSolidInterface_array_X, localSolidInterface_array_Y, localSolidInterface_array_Z = self.SolidSolver.getNodalInitialPositions()
            self.fillMatrix(localSolidInterface_array_X, localSolidInterface_array_Y, localSolidInterface_array_Z, 0)

        self.H.assemble()
        self.H_T.assemble()

        mpiPrint('\nInterpolation matrix is built.', mpiComm)

    def fillMatrix(self, solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc):
        """
        Des.
        """
            
        KDTree = self.createKDTree(solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z)
        localFluidInterface_array_X_init, localFluidInterface_array_Y_init, localFluidInterface_array_Z_init = self.FluidSolver.getNodalInitialPositions()
        for iVertex in range(self.nf_loc):
            posX = localFluidInterface_array_X_init[iVertex]
            posY = localFluidInterface_array_Y_init[iVertex]
            posZ = localFluidInterface_array_Z_init[iVertex]
            fluidPoint = np.array([posX, posY, posZ])
            distance, jVertex = KDTree.query(fluidPoint, 1)
            iGlobalVertexFluid = self.manager.getGlobalIndex('fluid', self.myid, iVertex)
            jGlobalVertexSolid = self.manager.getGlobalIndex('solid', iProc, jVertex)
            if distance > 1e-6:
                print("WARNING : Tolerance for matching meshes is not matched between node F{} and S{} : ({}, {}, {})<-->({}, {}, {}) , DISTANCE : {} !".format(iGlobalVertexFluid,jGlobalVertexSolid,posX, posY, posZ,solidInterfaceBuffRcv_X[jVertex], solidInterfaceBuffRcv_Y[jVertex], solidInterfaceBuffRcv_Z[jVertex], distance))
            self.H.setValue((iGlobalVertexFluid,jGlobalVertexSolid),1.0)
            self.H_T.setValue((jGlobalVertexSolid, iGlobalVertexFluid),1.0)

    def interpolateFluidLoadsOnSolidMesh(self):
        """
        Description
        """

        for iDim in range(3):
            self.H_T.mult(self.fluidInterfaceLoads.getData(iDim), self.solidInterfaceLoads.getData(iDim))

    def interpolateSolidDisplacementOnFluidMesh(self):
        """
        Description.
        """

        for iDim in range(3):
            self.H.mult(self.solidInterfaceDisplacement.getData(iDim), self.fluidInterfaceDisplacement.getData(iDim))

    def interpolateSolidHeatFluxOnFluidMesh(self):
        """
        Description.
        """

        for iDim in range(3):
            self.H.mult(self.solidInterfaceHeatFlux.getData(iDim), self.fluidInterfaceHeatFlux.getData(iDim))

    def interpolateSolidTemperatureOnFluidMesh(self):
        """
        Description
        """

        for iDim in range(1):
            self.H.mult(self.solidInterfaceTemperature.getData(iDim), self.fluidInterfaceTemperature.getData(iDim))

    def interpolateFluidHeatFluxOnSolidMesh(self):
        """
        Description.
        """

        for iDim in range(3):
            self.H_T.mult(self.fluidInterfaceHeatFlux.getData(iDim), self.solidInterfaceHeatFlux.getData(iDim))
        self.H_T.mult(self.fluidInterfaceNormalHeatFlux.getData(0), self.solidInterfaceNormalHeatFlux.getData(0))

    def interpolateFluidTemperatureOnSolidMesh(self):
        """
        Description.
        """

        for iDim in range(1):
            self.H_T.mult(self.fluidInterfaceTemperature.getData(iDim), self.solidInterfaceTemperature.getData(iDim))

    def interpolateFluidRobinTemperatureOnSolidMesh(self):
        """
        Des.
        """

        for iDim in range(1):
            self.H_T.mult(self.fluidInterfaceRobinTemperature.getData(iDim), self.solidInterfaceRobinTemperature.getData(iDim))
    

class RBFInterpolator(InterfaceInterpolator):
    """
    Description.
    """

    def __init__(self, Manager, FluidSolver, SolidSolver, RBFradius=0.1, mpiComm = None):
        """"
        Description.
        """

        InterfaceInterpolator.__init__(self, Manager, FluidSolver, SolidSolver, mpiComm)

        mpiPrint('\nSetting non matching meshes interpolator with Radial Basis Functions...', mpiComm)

        self.nf = self.manager.getNumberOfFluidInterfaceNodes()
        self.ns = self.manager.getNumberOfSolidInterfaceNodes()
        self.nf_loc = self.manager.getNumberOfLocalFluidInterfaceNodes()
        self.ns_loc = self.manager.getNumberOfLocalSolidInterfaceNodes()
        self.nDim = self.manager.getnDim()

        self.d = self.nDim+1
        self.radius = RBFradius

        self.generateInterfaceData()

        self.A = InterfaceMatrix((self.ns+self.d,self.ns+self.d), self.mpiComm)
        self.A_T = InterfaceMatrix((self.ns+self.d,self.ns+self.d), self.mpiComm)
        self.B = InterfaceMatrix((self.nf,self.ns+self.d), self.mpiComm)
        self.B_T = InterfaceMatrix((self.ns+self.d,self.nf), self.mpiComm)

        solidInterfaceProcessors = self.manager.getSolidInterfaceProcessors()
        fluidInterfaceProcessors = self.manager.getFluidInterfaceProcessors()
        solidPhysicalInterfaceNodesDistribution = self.manager.getSolidPhysicalInterfaceNodesDistribution()

        mpiPrint('\nBuilding interpolation matrices...', mpiComm)

        # Fill the matrix A
        if self.mpiComm != None:
            for iProc in solidInterfaceProcessors:
                if self.myid == iProc:
                    localSolidInterface_array_X, localSolidInterface_array_Y, localSolidInterface_array_Z = self.SolidSolver.getNodalInitialPositions()
                    for jProc in solidInterfaceProcessors:
                        self.mpiComm.Send(localSolidInterface_array_X, dest=jProc, tag=1)
                        self.mpiComm.Send(localSolidInterface_array_Y, dest=jProc, tag=2)
                        self.mpiComm.Send(localSolidInterface_array_Z, dest=jProc, tag=3)
                if self.myid in solidInterfaceProcessors:
                    sizeOfBuff = solidPhysicalInterfaceNodesDistribution[iProc]
                    solidInterfaceBuffRcv_X = np.zeros(sizeOfBuff)
                    solidInterfaceBuffRcv_Y = np.zeros(sizeOfBuff)
                    solidInterfaceBuffRcv_Z = np.zeros(sizeOfBuff)
                    self.mpiComm.Recv(solidInterfaceBuffRcv_X, iProc, tag=1)
                    self.mpiComm.Recv(solidInterfaceBuffRcv_Y, iProc, tag=2)
                    self.mpiComm.Recv(solidInterfaceBuffRcv_Z, iProc, tag=3)
                    self.fillMatrixA(solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc)
        else:
            localSolidInterface_array_X, localSolidInterface_array_Y, localSolidInterface_array_Z = self.SolidSolver.getNodalInitialPositions()
            self.fillMatrixA(localSolidInterface_array_X, localSolidInterface_array_Y, localSolidInterface_array_Z, 0)

        self.A.assemble()
        self.A_T.assemble()

        mpiPrint('\nMatrix A is built.', mpiComm)

        # Fill the matrix B
        if self.mpiComm != None:
            for iProc in solidInterfaceProcessors:
                if self.myid == iProc:
                    for jProc in fluidInterfaceProcessors:
                        self.mpiComm.Send(localSolidInterface_array_X, dest=jProc, tag=1)
                        self.mpiComm.Send(localSolidInterface_array_Y, dest=jProc, tag=2)
                        self.mpiComm.Send(localSolidInterface_array_Z, dest=jProc, tag=3)
                if self.myid in fluidInterfaceProcessors:
                    sizeOfBuff = solidPhysicalInterfaceNodesDistribution[iProc]
                    solidInterfaceBuffRcv_X = np.zeros(sizeOfBuff)
                    solidInterfaceBuffRcv_Y = np.zeros(sizeOfBuff)
                    solidInterfaceBuffRcv_Z = np.zeros(sizeOfBuff)
                    self.mpiComm.Recv(solidInterfaceBuffRcv_X, iProc, tag=1)
                    self.mpiComm.Recv(solidInterfaceBuffRcv_Y, iProc, tag=2)
                    self.mpiComm.Recv(solidInterfaceBuffRcv_Z, iProc, tag=3)
                    self.fillMatrixB(solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc)
        else:
            self.fillMatrixB(localSolidInterface_array_X, localSolidInterface_array_Y, localSolidInterface_array_Z, 0)

        self.B.assemble()
        self.B_T.assemble()

        mpiPrint('\nMatrix B is built.', mpiComm)

    def fillMatrixA(self, solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc):
        """
        Description.
        """

        KDTree = self.createKDTree(solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z)
        localSolidInterface_array_X_init, localSolidInterface_array_Y_init, localSolidInterface_array_Z_init = self.SolidSolver.getNodalInitialPositions()
        for iVertex in range(self.ns_loc):
            posX = localSolidInterface_array_X_init[iVertex]
            posY = localSolidInterface_array_Y_init[iVertex]
            posZ = localSolidInterface_array_Z_init[iVertex]
            solidPoint = np.array([posX, posY, posZ])
            solidVertices = KDTree.query_ball_point(solidPoint, 1.01*self.radius)
            iGlobalVertexSolid = self.manager.getGlobalIndex('solid', self.myid, iVertex)
            for jVertex in solidVertices:
                jGlobalVertexSolid = self.manager.getGlobalIndex('solid', iProc, jVertex)
                solidQuery = np.array([solidInterfaceBuffRcv_X[jVertex], solidInterfaceBuffRcv_Y[jVertex], solidInterfaceBuffRcv_Z[jVertex]])
                distance = self.distance(solidPoint, solidQuery)
                phi = self.__PHI(distance)
                self.A.setValue((iGlobalVertexSolid, jGlobalVertexSolid), phi)
                self.A_T.setValue((jGlobalVertexSolid, iGlobalVertexSolid), phi)
            self.A.setValue((iGlobalVertexSolid, self.ns), 1.0)
            self.A.setValue((iGlobalVertexSolid, self.ns+1), posX)
            self.A.setValue((iGlobalVertexSolid, self.ns+2), posY)
            self.A_T.setValue((self.ns, iGlobalVertexSolid), 1.0)
            self.A_T.setValue((self.ns+1, iGlobalVertexSolid), posX)
            self.A_T.setValue((self.ns+2, iGlobalVertexSolid), posY)
            if self.nDim == 3:
                self.A.setValue((iGlobalVertexSolid, self.ns+3), posZ)
                self.A_T.setValue((self.ns+3, iGlobalVertexSolid), posZ)

    def fillMatrixB(self, solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc):
        """
        Description.
        """

        KDTree = self.createKDTree(solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z)
        localFluidInterface_array_X_init, localFluidInterface_array_Y_init, localFluidInterface_array_Z_init = self.FluidSolver.getNodalInitialPositions()
        for iVertex in range(self.nf_loc):
            posX = localFluidInterface_array_X_init[iVertex]
            posY = localFluidInterface_array_Y_init[iVertex]
            posZ = localFluidInterface_array_Z_init[iVertex]
            fluidPoint = np.array([posX, posY, posZ])
            solidVertices = KDTree.query_ball_point(fluidPoint, 1.01*self.radius)
            iGlobalVertexFluid = self.manager.getGlobalIndex('fluid', self.myid, iVertex)
            for jVertex in solidVertices:
                jGlobalVertexSolid = self.manager.getGlobalIndex('solid', iProc, jVertex)
                solidQuery = np.array([solidInterfaceBuffRcv_X[jVertex], solidInterfaceBuffRcv_Y[jVertex], solidInterfaceBuffRcv_Z[jVertex]])
                distance = self.distance(fluidPoint, solidQuery)
                phi = self.__PHI(distance)
                self.B.setValue((iGlobalVertexFluid, jGlobalVertexSolid), phi)
                self.B_T.setValue((jGlobalVertexSolid, iGlobalVertexFluid), phi)
            self.B.setValue((iGlobalVertexFluid, self.ns), 1.0)
            self.B.setValue((iGlobalVertexFluid, self.ns+1), posX)
            self.B.setValue((iGlobalVertexFluid, self.ns+2), posY)
            self.B_T.setValue((self.ns, iGlobalVertexFluid), 1.0)
            self.B_T.setValue((self.ns+1, iGlobalVertexFluid), posX)
            self.B_T.setValue((self.ns+2, iGlobalVertexFluid), posY)
            if self.nDim == 3:
                self.B.setValue((iGlobalVertexFluid, self.ns+3), posZ)
                self.B_T.setValue((self.ns+3, iGlobalVertexFluid), posZ)

    def __PHI(self, distance):
        """
        Description.
        """

        phi = 0.0
        eps = distance/self.radius

        if eps < 1:
            phi = ((1.0-eps)**4)*(4.0*eps+1.0)  #BW - CPC4
            #phi = ((1.0-eps)**2)*(0.5*eps+1.0)   #EU - CPC2
        else:
            phi = 0.0

        return phi

    def interpolateFluidLoadsOnSolidMesh(self):
        """
        Description
        """

        gamma_array = FlexInterfaceData(self.ns + self.d, 3, self.mpiComm)
        Solver = LinearSolver(self.A_T, self.mpiComm)

        for iDim in range(3):
            self.B_T.mult(self.fluidInterfaceLoads.getData(iDim), gamma_array.getData(iDim))

        for iDim in range(3):
            Solver.solve(gamma_array.getData(iDim), self.solidInterfaceLoads.getData(iDim))

    def interpolateSolidDisplacementOnFluidMesh(self):
        """
        Description.
        """

        gamma_array = FlexInterfaceData(self.ns + self.d, 3, self.mpiComm)
        Solver = LinearSolver(self.A, self.mpiComm)

        for iDim in range(3):
            Solver.solve(self.solidInterfaceDisplacement.getData(iDim), gamma_array.getData(iDim))

        for iDim in range(3):
            self.B.mult(gamma_array.getData(iDim), self.fluidInterfaceDisplacement.getData(iDim))

class TPSInterpolator(RBFInterpolator):
    """
    Des.
    """

    def __init__(self, Manager, FluidSolver, SolidSolver, mpiComm=None):
        """
        des.
        """

        InterfaceInterpolator.__init__(self, Manager, FluidSolver, SolidSolver, mpiComm)

        mpiPrint('\nSetting non matching meshes interpolator with Thin Plate Spline...', mpiComm)

        self.nf = self.manager.getNumberOfFluidInterfaceNodes()
        self.ns = self.manager.getNumberOfSolidInterfaceNodes()
        self.nf_loc = self.manager.getNumberOfLocalFluidInterfaceNodes()
        self.ns_loc = self.manager.getNumberOfLocalSolidInterfaceNodes()
        self.nDim = self.manager.getnDim()

        self.d = self.nDim+1

        self.generateInterfaceData()

        self.A = InterfaceMatrix((self.ns+self.d,self.ns+self.d), self.mpiComm)
        self.A_T = InterfaceMatrix((self.ns+self.d,self.ns+self.d), self.mpiComm)
        self.B = InterfaceMatrix((self.nf,self.ns+self.d), self.mpiComm)
        self.B_T = InterfaceMatrix((self.ns+self.d,self.nf), self.mpiComm)

        solidInterfaceProcessors = self.manager.getSolidInterfaceProcessors()
        fluidInterfaceProcessors = self.manager.getFluidInterfaceProcessors()
        solidPhysicalInterfaceNodesDistribution = self.manager.getSolidPhysicalInterfaceNodesDistribution()

        mpiPrint('\nBuilding interpolation matrices...', mpiComm)

        # Fill the matrix A
        if self.mpiComm != None:
            for iProc in solidInterfaceProcessors:
                if self.myid == iProc:
                    localSolidInterface_array_X, localSolidInterface_array_Y, localSolidInterface_array_Z = self.SolidSolver.getNodalInitialPositions()
                    for jProc in solidInterfaceProcessors:
                        self.mpiComm.Send(localSolidInterface_array_X, dest=jProc, tag=1)
                        self.mpiComm.Send(localSolidInterface_array_Y, dest=jProc, tag=2)
                        self.mpiComm.Send(localSolidInterface_array_Z, dest=jProc, tag=3)
                if self.myid in solidInterfaceProcessors:
                    sizeOfBuff = solidPhysicalInterfaceNodesDistribution[iProc]
                    solidInterfaceBuffRcv_X = np.zeros(sizeOfBuff)
                    solidInterfaceBuffRcv_Y = np.zeros(sizeOfBuff)
                    solidInterfaceBuffRcv_Z = np.zeros(sizeOfBuff)
                    self.mpiComm.Recv(solidInterfaceBuffRcv_X, iProc, tag=1)
                    self.mpiComm.Recv(solidInterfaceBuffRcv_Y, iProc, tag=2)
                    self.mpiComm.Recv(solidInterfaceBuffRcv_Z, iProc, tag=3)
                    self.fillMatrixA(solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc)
        else:
            localSolidInterface_array_X, localSolidInterface_array_Y, localSolidInterface_array_Z = self.SolidSolver.getNodalInitialPositions()
            self.fillMatrixA(localSolidInterface_array_X, localSolidInterface_array_Y, localSolidInterface_array_Z, 0)

        self.A.assemble()
        self.A_T.assemble()

        mpiPrint('\nMatrix A is built.', mpiComm)

        # Fill the matrix B
        if self.mpiComm != None:
            for iProc in solidInterfaceProcessors:
                if self.myid == iProc:
                    for jProc in fluidInterfaceProcessors:
                        self.mpiComm.Send(localSolidInterface_array_X, dest=jProc, tag=1)
                        self.mpiComm.Send(localSolidInterface_array_Y, dest=jProc, tag=2)
                        self.mpiComm.Send(localSolidInterface_array_Z, dest=jProc, tag=3)
                if self.myid in fluidInterfaceProcessors:
                    sizeOfBuff = solidPhysicalInterfaceNodesDistribution[iProc]
                    solidInterfaceBuffRcv_X = np.zeros(sizeOfBuff)
                    solidInterfaceBuffRcv_Y = np.zeros(sizeOfBuff)
                    solidInterfaceBuffRcv_Z = np.zeros(sizeOfBuff)
                    self.mpiComm.Recv(solidInterfaceBuffRcv_X, iProc, tag=1)
                    self.mpiComm.Recv(solidInterfaceBuffRcv_Y, iProc, tag=2)
                    self.mpiComm.Recv(solidInterfaceBuffRcv_Z, iProc, tag=3)
                    self.fillMatrixB(solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc)
        else:
            self.fillMatrixB(localSolidInterface_array_X, localSolidInterface_array_Y, localSolidInterface_array_Z, 0)

        self.B.assemble()
        self.B_T.assemble()

        mpiPrint('\nMatrix B is built.', mpiComm)

    def fillMatrixA(self, solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc):
        """
        Description.
        """

        localSolidInterface_array_X_init, localSolidInterface_array_Y_init, localSolidInterface_array_Z_init = self.SolidSolver.getNodalInitialPositions()
        for iVertex in range(self.ns_loc):
            posX = localSolidInterface_array_X_init[iVertex]
            posY = localSolidInterface_array_Y_init[iVertex]
            posZ = localSolidInterface_array_Z_init[iVertex]
            solidPoint = np.array([posX, posY, posZ])
            iGlobalVertexSolid = self.manager.getGlobalIndex('solid', self.myid, iVertex)
            for jVertex in range(solidInterfaceBuffRcv_X.shape[0]):
                jGlobalVertexSolid = self.manager.getGlobalIndex('solid', iProc, jVertex)
                solidQuery = np.array([solidInterfaceBuffRcv_X[jVertex], solidInterfaceBuffRcv_Y[jVertex], solidInterfaceBuffRcv_Z[jVertex]])
                distance = self.distance(solidPoint, solidQuery)
                phi = self.__PHI(distance)
                self.A.setValue((iGlobalVertexSolid, jGlobalVertexSolid), phi)
                #print("Set : {}<-->{} = {}".format(iGlobalVertexSolid, jGlobalVertexSolid,phi))
                self.A_T.setValue((jGlobalVertexSolid, iGlobalVertexSolid), phi)
            self.A.setValue((iGlobalVertexSolid, self.ns), 1.0)
            self.A.setValue((iGlobalVertexSolid, self.ns+1), posX)
            self.A.setValue((iGlobalVertexSolid, self.ns+2), posY)
            self.A_T.setValue((self.ns, iGlobalVertexSolid), 1.0)
            self.A_T.setValue((self.ns+1, iGlobalVertexSolid), posX)
            self.A_T.setValue((self.ns+2, iGlobalVertexSolid), posY)
            if self.nDim == 3:
                self.A.setValue((iGlobalVertexSolid, self.ns+3), posZ)
                self.A_T.setValue((self.ns+3, iGlobalVertexSolid), posZ)

    def fillMatrixB(self, solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc):
        """
        Description.
        """

        localFluidInterface_array_X_init, localFluidInterface_array_Y_init, localFluidInterface_array_Z_init = self.FluidSolver.getNodalInitialPositions()
        for iVertex in range(self.nf_loc):
            posX = localFluidInterface_array_X_init[iVertex]
            posY = localFluidInterface_array_Y_init[iVertex]
            posZ = localFluidInterface_array_Z_init[iVertex]
            fluidPoint = np.array([posX, posY, posZ])
            iGlobalVertexFluid = self.manager.getGlobalIndex('fluid', self.myid, iVertex)
            for jVertex in range(solidInterfaceBuffRcv_X.shape[0]):
                jGlobalVertexSolid = self.manager.getGlobalIndex('solid', iProc, jVertex)
                solidQuery = np.array([solidInterfaceBuffRcv_X[jVertex], solidInterfaceBuffRcv_Y[jVertex], solidInterfaceBuffRcv_Z[jVertex]])
                distance = self.distance(fluidPoint, solidQuery)
                phi = self.__PHI(distance)
                self.B.setValue((iGlobalVertexFluid, jGlobalVertexSolid), phi)
                self.B_T.setValue((jGlobalVertexSolid, iGlobalVertexFluid), phi)
            self.B.setValue((iGlobalVertexFluid, self.ns), 1.0)
            self.B.setValue((iGlobalVertexFluid, self.ns+1), posX)
            self.B.setValue((iGlobalVertexFluid, self.ns+2), posY)
            self.B_T.setValue((self.ns, iGlobalVertexFluid), 1.0)
            self.B_T.setValue((self.ns+1, iGlobalVertexFluid), posX)
            self.B_T.setValue((self.ns+2, iGlobalVertexFluid), posY)
            if self.nDim == 3:
                self.B.setValue((iGlobalVertexFluid, self.ns+3), posZ)
                self.B_T.setValue((self.ns+3, iGlobalVertexFluid), posZ)

    def __PHI(self, distance):
        """
        Description.
        """

        phi = 0.0

        if distance > 0.0:
            phi = (distance**2)*np.log10(distance)
        else:
            phi = 0.0

        return phi

# ----------------------------------------------------------------------
#    Algorithm class
# ----------------------------------------------------------------------

class Algortihm:
    """
    Des.
    """
    
    def __init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, Criterion, nbFSIIterMax, deltaT, totTime, timeIterTreshold=-1, mpiComm=None):
        """
        Des.
        """
        
        mpiPrint('\n***************************** Initializing FSI algorithm *****************************', mpiComm)
        
        self.run_t0 = 0. # timer to estimate run time: initial time
        self.run_tf = 0. # timer to estimate run time: final time
        
        self.mpiComm = mpiComm
        self.manager = Manager
        self.FluidSolver = FluidSolver
        self.SolidSolver = SolidSolver
        self.interfaceInterpolator = InterfaceInterpolator
        self.criterion = Criterion

        self.nbFSIIterMax = nbFSIIterMax        
        self.deltaT = deltaT
        self.totTime = totTime
        self.timeIterTreshold = timeIterTreshold
        self.writeInFSIloop = False
        
        self.FSIIter = 0
        self.errValue = 0.0
        self.FSIConv = False
        self.totNbOfFSIIt = 0
        
        if self.mpiComm != None:
            self.myid = self.mpiComm.Get_rank()
            self.mpiSize = self.mpiComm.Get_size()
        else:
            self.myid = 0
            self.mpiSize = 1
        
        self.alpha_0 = 1.0
        self.alpha_1 = 0.5
        
        ns = self.interfaceInterpolator.getNs()
        d = self.interfaceInterpolator.getd()
        
        self.solidInterfaceVelocity = FlexInterfaceData(ns+d, 3, self.mpiComm)
        self.solidInterfaceVelocitynM1 = FlexInterfaceData(ns+d, 3, self.mpiComm)

        self.solidInterfaceResidual = FlexInterfaceData(ns+d, 3, self.mpiComm)

        self.solidHeatFluxResidual = None
        self.solidTemperatureResidual = None
        if self.manager.withCht:
            self.solidHeatFluxResidual = FlexInterfaceData(ns+d, 3, self.mpiComm)
            self.solidTemperatureResidual = FlexInterfaceData(ns+d, 1, self.mpiComm)
    
    def setFSIInitialConditions(self, time):
        """
        Des.
        """
        
        if self.manager.computationType == 'unsteady':
            if self.myid in self.manager.getSolidSolverProcessors():
                self.SolidSolver.setInitialDisplacements()
            self.interfaceInterpolator.getDisplacementFromSolidSolver()
            self.interfaceInterpolator.interpolateSolidDisplacementOnFluidMesh()
            self.interfaceInterpolator.setDisplacementToFluidSolver(time)
            self.FluidSolver.setInitialMeshDeformation()
            if self.manager.withCht:
                if self.interfaceInterpolator.chtTransferMethod == 'hFFB' or self.interfaceInterpolator.chtTransferMethod == 'TFFB':
                    #self.interfaceInterpolator.getHeatFluxFromSolidSolver()
                    #self.interfaceInterpolator.interpolateSolidHeatFluxOnFluidMesh()
                    #self.interfaceInterpolator.setHeatFluxToFluidSolver(time)
                    self.FluidSolver.setInitialInterfaceHeatFlux()
                elif self.interfaceInterpolator.chtTransferMethod == 'hFTB' or self.interfaceInterpolator.chtTransferMethod == 'FFTB':
                    self.FluidSolver.setInitialInterfaceTemperature()
        else:
            self.interfaceInterpolator.getDisplacementFromSolidSolver()
            if self.manager.withCht:
                if self.interfaceInterpolator.chtTransferMethod == 'hFFB' or self.interfaceInterpolator.chtTransferMethod == 'TFFB':
                    #self.interfaceInterpolator.getHeatFluxFromSolidSolver()
                    #self.interfaceInterpolator.interpolateSolidHeatFluxOnFluidMesh()
                    #self.interfaceInterpolator.setHeatFluxToFluidSolver(time)
                    self.FluidSolver.setInitialInterfaceHeatFlux()
                elif self.interfaceInterpolator.chtTransferMethod == 'hFTB' or self.interfaceInterpolator.chtTransferMethod == 'FFTB':
                    self.FluidSolver.setInitialInterfaceTemperature()
    
    def computeSolidInterfaceResidual(self):
        """
        Des.
        """
        
        ns = self.interfaceInterpolator.getNs()
        d = self.interfaceInterpolator.getd()
        
        # --- Get the predicted (computed) solid interface displacement from the solid solver --- #
        predictedDisplacement = FlexInterfaceData(ns+d, 3, self.mpiComm)
        
        if self.myid in self.manager.getSolidInterfaceProcessors():
            localSolidInterfaceDisp_X, localSolidInterfaceDisp_Y, localSolidInterfaceDisp_Z = self.SolidSolver.getNodalDisplacements()
            for iVertex in range(self.manager.getNumberOfLocalSolidInterfaceNodes()):
                iGlobalVertex = self.manager.getGlobalIndex('solid', self.myid, iVertex)
                predictedDisplacement[iGlobalVertex] = [localSolidInterfaceDisp_X[iVertex], localSolidInterfaceDisp_Y[iVertex], localSolidInterfaceDisp_Z[iVertex]]
        
        predictedDisplacement.assemble()
        
        # --- Calculate the residual (vector and norm) --- #
        mpiPrint("\nCompute FSI residual based on solid interface displacement.", self.mpiComm)
        self.solidInterfaceResidual = predictedDisplacement - self.interfaceInterpolator.solidInterfaceDisplacement
        
        return self.solidInterfaceResidual

    def computeSolidInterfaceResidual_CHT(self):
        """
        Des.
        """

        ns = self.interfaceInterpolator.getNs()
        d = self.interfaceInterpolator.getd()

        predictedHF = FlexInterfaceData(ns+d, 3, self.mpiComm)
        predictedTemp = FlexInterfaceData(ns+d, 1, self.mpiComm)

        if self.myid in self.manager.getSolidInterfaceProcessors():
            localSolidInterfaceHeatFlux_X, localSolidInterfaceHeatFlux_Y, localSolidInterfaceHeatFlux_Z = self.SolidSolver.getNodalHeatFluxes()
            localSolidInterfaceTemperature = self.SolidSolver.getNodalTemperatures()
            for iVertex in range(self.manager.getNumberOfLocalSolidInterfaceNodes()):
                iGlobalVertex = self.manager.getGlobalIndex('solid', self.myid, iVertex)
                predictedHF[iGlobalVertex] = [localSolidInterfaceHeatFlux_X[iVertex], localSolidInterfaceHeatFlux_Y[iVertex], localSolidInterfaceHeatFlux_Z[iVertex]]
                predictedTemp[iGlobalVertex] = [localSolidInterfaceTemperature[iVertex]]

        predictedHF.assemble()
        predictedTemp.assemble()

        self.solidHeatFluxResidual = predictedHF - self.interfaceInterpolator.solidInterfaceHeatFlux
        self.solidTemperatureResidual = predictedTemp - self.interfaceInterpolator.solidInterfaceTemperature

        #return self.solidHeatFluxResidual

        if self.interfaceInterpolator.chtTransferMethod == 'hFFB' or self.interfaceInterpolator.chtTransferMethod == 'TFFB':
            mpiPrint("\nCompute CHT residual based on solid interface heat flux.", self.mpiComm)
            return self.solidHeatFluxResidual
        elif self.interfaceInterpolator.chtTransferMethod == 'hFTB' or self.interfaceInterpolator.chtTransferMethod == 'FFTB':
            mpiPrint("\nCompute CHT residual based on solid interface temperature.", self.mpiComm)
            return self.solidTemperatureResidual
        else:
            return None
    
    def solidDisplacementPredictor(self):
        """
        Des
        """
        
        # --- Get the velocity (current and previous time step) of the solid interface from the solid solver --- #
        if self.myid in self.manager.getSolidInterfaceProcessors():
            localSolidInterfaceVel_X, localSolidInterfaceVel_Y, localSolidInterfaceVel_Z = self.SolidSolver.getNodalVelocity()
            localSolidInterfaceVelNm1_X, localSolidInterfaceVelNm1_Y, localSolidInterfaceVelNm1_Z = self.SolidSolver.getNodalVelocityNm1()
            for iVertex in range(self.manager.getNumberOfLocalSolidInterfaceNodes()):
                iGlobalVertex = self.manager.getGlobalIndex('solid', self.myid, iVertex)
                self.solidInterfaceVelocity[iGlobalVertex] = [localSolidInterfaceVel_X[iVertex], localSolidInterfaceVel_Y[iVertex], localSolidInterfaceVel_Z[iVertex]]
                self.solidInterfaceVelocitynM1[iGlobalVertex] = [localSolidInterfaceVelNm1_X[iVertex], localSolidInterfaceVelNm1_Y[iVertex], localSolidInterfaceVelNm1_Z[iVertex]]
        
        self.solidInterfaceVelocity.assemble()
        self.solidInterfaceVelocitynM1.assemble()
        
        # --- Predict the solid position for the next time step --- #
        self.interfaceInterpolator.solidInterfaceDisplacement += (self.alpha_0*self.deltaT*self.solidInterfaceVelocity + self.alpha_1*self.deltaT*(self.solidInterfaceVelocity-self.solidInterfaceVelocitynM1))

    def iniRealTimeData(self):
        """
        Des
        """
        
        if self.myid in self.manager.getSolidSolverProcessors():
            self.SolidSolver.initRealTimeData()
        histFile = open('FSIhistory.ascii', "w")
        histFile.write("TimeIter\tTime\tFSIError\tCHTError\tFSINbIter\n")
        histFile.close()
    
    def writeRealTimeData(self, timeIter, time, error, chtError=0.0):
        """
        Des
        """
        
        if self.myid == 0:
            self.FluidSolver.saveRealTimeData(time, self.FSIIter)
            if timeIter >= self.timeIterTreshold:
                self.SolidSolver.saveRealTimeData(time, self.FSIIter)
            histFile = open('FSIhistory.ascii', "a")
            histFile.write(str(timeIter) + '\t' + str(time) + '\t' + str(error) + '\t' + str(chtError) + '\t' + str(self.FSIIter) + '\n')
            histFile.close()
    
    def getMeanNbOfFSIIt(self, timeIter):
        """
        Des
        """
        
        return float(self.totNbOfFSIIt)/timeIter
    
    def printExitInfo(self, timeIter, time, error):
        """
        Des
        """
        
        mpiPrint('[cpu FSI]: ' + str(self.run_tf - self.run_t0) + 's', self.mpiComm)
        mpiPrint('[Time steps FSI]: ' + str(timeIter), self.mpiComm)
        mpiPrint('[Successful Run FSI]: ' + str(time >= self.totTime - self.deltaT), self.mpiComm)
        mpiPrint('[Mean n. of FSI Iterations]: ' + str(self.getMeanNbOfFSIIt(timeIter)), self.mpiComm)
        
        if self.myid == 0 :
            self.FluidSolver.printRealTimeData(time, self.FSIIter)
            self.SolidSolver.printRealTimeData(time, self.FSIIter)
        
        mpiPrint('RES-FSI-FSIhistory: ' + str(timeIter) + '\t' + str(time) + '\t' + str(error) + '\t' + str(self.FSIIter) + '\n', self.mpiComm)
    
    def run(self):
        """
        Des.
        """
        
        # --- Initialize output manager --- #
        self.iniRealTimeData()
        
        mpiPrint('\n**********************************', self.mpiComm)
        mpiPrint('*         Begin FSI computation            *', self.mpiComm)
        mpiPrint('**********************************\n', self.mpiComm)
        
        self.run_t0 = tm.time()
        
        #If no restart
        mpiPrint('Setting FSI initial conditions...', self.mpiComm)
        self.setFSIInitialConditions(0.0)
        mpiPrint('\nFSI initial conditions are set', self.mpiComm)
        
        if self.manager.computationType == 'unsteady':
            self.__unsteadyRun()
        else:
            time = self.totTime
            timeIter = 1
            self.deltaT = self.totTime
            self.writeInFSIloop = True
            self.fsiCoupling(timeIter, time)
            self.totNbOfFSIIt = self.FSIIter
            self.run_tf = tm.time()
            self.printExitInfo(timeIter, time, self.errValue)
            
        mpiBarrier(self.mpiComm)
        
        mpiPrint('\n*************************', self.mpiComm)
        mpiPrint('*    End FSI computation    *', self.mpiComm)
        mpiPrint('*************************\n', self.mpiComm)
    
    def __unsteadyRun(self):
        """
        Des.
        """
        
        #If no restart
        nbTimeIter = int((self.totTime/self.deltaT)-1)
        time = 0.0
        timeIter = 0
        
        mpiPrint('Begin time integration\n', self.mpiComm)
        
        # --- External temporal loop --- #
        while timeIter <= nbTimeIter:
            
            mpiPrint("\n>>>> Time iteration {} <<<<".format(timeIter), self.mpiComm)
            
            # --- Preprocess the temporal iteration --- #
            self.FluidSolver.preprocessTimeIter(timeIter)
            if self.myid in self.manager.getSolidSolverProcessors():
                self.SolidSolver.preprocessTimeIter(timeIter)
            
            # --- Internal FSI loop --- #
            self.fsiCoupling(timeIter, time)
            # --- End of FSI loop --- #
            
            mpiBarrier(self.mpiComm)
            
            self.totNbOfFSIIt += self.FSIIter
            
            # --- Update the fluid and solid solver for the next time step --- #
            if self.myid in self.manager.getSolidSolverProcessors():
                self.SolidSolver.update()
            self.FluidSolver.update(self.deltaT)
            
            # --- Write fluid and solid solution, update FSI history  ---#
            self.FluidSolver.save(timeIter)
            
            if self.myid in self.manager.getSolidSolverProcessors():
                self.SolidSolver.save()
            
            self.writeRealTimeData(timeIter, time, self.errValue, self.errValue_CHT)
            
            # --- Perform some remeshing if necessary
            if self.myid in self.manager.getSolidSolverProcessors():
                self.SolidSolver.remeshing()
            self.FluidSolver.remeshing()
            # ---
            
            if timeIter >= self.timeIterTreshold:
                # --- Displacement predictor for the next time step and update of the solid solution --- #
                mpiPrint('\nSolid displacement prediction for next time step', self.mpiComm)
                self.solidDisplacementPredictor()
            
            timeIter += 1
            time += self.deltaT
        # --- End of the temporal loop --- #
        
        self.run_tf = tm.time()
        self.printExitInfo(timeIter, time, self.errValue)
    
class AlgortihmBGSStaticRelax(Algortihm):
    """
    Des.
    """
    
    def __init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, Criterion, nbFSIIterMax, deltaT, totTime, timeIterTreshold=-1, omegaMax=1.0, mpiComm=None):
        """
        Des.
        """
        
        Algortihm.__init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, Criterion, nbFSIIterMax, deltaT, totTime, timeIterTreshold, mpiComm)
        
        self.omegaMax = omegaMax
        self.omega = omegaMax
    
    def fsiCoupling(self, timeIter, time):
        """
        Block Gauss Seidel (BGS) method for strong coupling FSI
        """
        
        if timeIter > self.timeIterTreshold:
            nbFSIIter = self.nbFSIIterMax
            mpiPrint('\n*************** Enter Block Gauss Seidel (BGS) method for strong coupling FSI ***************', self.mpiComm)
        else:
             nbFSIIter = 1
        
        self.FSIIter = 0
        self.FSIConv = False
        self.errValue = 1.0
        self.errValue_CHT = 1.0

        solidHasRun = False
        
        while ((self.FSIIter < nbFSIIter) and (not self.criterion.isVerified(self.errValue, self.errValue_CHT))):
            mpiPrint("\n>>>> FSI iteration {} <<<<\n".format(self.FSIIter), self.mpiComm)
            
            # --- Mesh morphing step (displacements interpolation, displacements communication, and mesh morpher call) --- #
            self.interfaceInterpolator.interpolateSolidDisplacementOnFluidMesh()
            self.interfaceInterpolator.setDisplacementToFluidSolver(time)
            mpiPrint('\nPerforming mesh deformation...\n', self.mpiComm)
            self.FluidSolver.meshUpdate(timeIter)
            # --- CHT data transfer (heat flux or temperature) --- 
            if solidHasRun:
                if self.interfaceInterpolator.chtTransferMethod == 'TFFB' or self.interfaceInterpolator.chtTransferMethod == 'hFFB':
                    self.interfaceInterpolator.interpolateSolidHeatFluxOnFluidMesh()
                    self.interfaceInterpolator.setHeatFluxToFluidSolver(time)
                elif self.interfaceInterpolator.chtTransferMethod == 'FFTB' or self.interfaceInterpolator.chtTransferMethod == 'hFTB':
                    self.interfaceInterpolator.interpolateSolidTemperatureOnFluidMesh()
                    self.interfaceInterpolator.setTemperatureToFluidSolver(time)
            
            # --- Fluid solver call for FSI subiteration --- #
            mpiPrint('\nLaunching fluid solver...', self.mpiComm)
            self.FluidSolver.run(time-self.deltaT, time)
            mpiBarrier(self.mpiComm)
            
            # --- Surface fluid loads interpolation and communication --- #
            mpiPrint('\nProcessing interface fluid loads...\n', self.mpiComm)
            self.interfaceInterpolator.getLoadsFromFluidSolver()
            if self.interfaceInterpolator.chtTransferMethod == 'TFFB':
                self.interfaceInterpolator.getTemperatureFromFluidSolver()
            elif self.interfaceInterpolator.chtTransferMethod == 'FFTB':
                self.interfaceInterpolator.getHeatFluxFromFluidSolver()
            elif self.interfaceInterpolator.chtTransferMethod == 'hFFB' or self.interfaceInterpolator.chtTransferMethod == 'hFTB':
                self.interfaceInterpolator.getRobinTemperatureFromFluidSolver()
            mpiBarrier(self.mpiComm)
            if timeIter > self.timeIterTreshold:
                self.interfaceInterpolator.interpolateFluidLoadsOnSolidMesh()
                self.interfaceInterpolator.setLoadsToSolidSolver(time)
                if self.interfaceInterpolator.chtTransferMethod == 'TFFB':
                    self.interfaceInterpolator.interpolateFluidTemperatureOnSolidMesh()
                    self.interfaceInterpolator.setTemperatureToSolidSolver(time)
                elif self.interfaceInterpolator.chtTransferMethod == 'hFFB' or self.interfaceInterpolator.chtTransferMethod == 'hFTB':
                    self.interfaceInterpolator.interpolateFluidRobinTemperatureOnSolidMesh()
                    self.interfaceInterpolator.setRobinHeatFluxToSolidSolver(time)
                elif self.interfaceInterpolator.chtTransferMethod == 'FFTB':
                    self.interfaceInterpolator.interpolateFluidHeatFluxOnSolidMesh()
                    self.interfaceInterpolator.setHeatFluxToSolidSolver(time)
                
                # --- Solid solver call for FSI subiteration --- #
                mpiPrint('\nLaunching solid solver...\n', self.mpiComm)
                if self.myid in self.manager.getSolidSolverProcessors():
                    self.SolidSolver.run(time-self.deltaT, time)
                solidHasRun = True

                # --- Compute and monitor the FSI residual --- #
                res = self.computeSolidInterfaceResidual()
                res_CHT = self.computeSolidInterfaceResidual_CHT()
                self.errValue = self.criterion.update(res)
                self.errValue_CHT = self.criterion.updateHeatFlux(res_CHT)
                mpiPrint('\nFSI error value : {}\n'.format(self.errValue), self.mpiComm)
                mpiPrint('\nCHT error value : {}\n'.format(self.errValue_CHT), self.mpiComm)
                self.FSIConv = self.criterion.isVerified(self.errValue, self.errValue_CHT)
                
                # --- Relaxe the solid position --- #
                mpiPrint('\nProcessing interface displacements...\n', self.mpiComm)
                self.relaxSolidPosition()
                self.relaxCHT()
            
            if self.writeInFSIloop == True:
                self.writeRealTimeData(timeIter, time, self.errValue, self.errValue_CHT)
            
            self.FSIIter += 1
            if self.manager.computationType != 'unsteady':
                time += self.deltaT

            # --- Update the solvers for the next BGS iteration --- #
            if self.myid in self.manager.getSolidSolverProcessors():
                self.SolidSolver.bgsUpdate()
            self.FluidSolver.bgsUpdate()
        
        if timeIter > self.timeIterTreshold:
            mpiPrint('\n*************** BGS is converged ***************', self.mpiComm)
    
    def setOmega(self):
        """
        Des.
        """
        
        self.omega = self.omegaMax
        
        mpiPrint('Static under-relaxation step with parameter {}'.format(self.omega), self.mpiComm)
    
    def relaxSolidPosition(self):
        """
        Des.
        """
        
        # --- Set the relaxation parameter --- #
        self.setOmega()
        
        # --- Relax the solid interface position --- #
        self.interfaceInterpolator.solidInterfaceDisplacement += (self.omega*self.solidInterfaceResidual)

    def relaxCHT(self):
        """
        Des.
        """

        if self.interfaceInterpolator.chtTransferMethod == 'hFFB' or self.interfaceInterpolator.chtTransferMethod == 'TFFB':
            self.interfaceInterpolator.solidInterfaceHeatFlux += self.solidHeatFluxResidual
        elif self.interfaceInterpolator.chtTransferMethod == 'hFTB' or self.interfaceInterpolator.chtTransferMethod == 'FFTB':
            self.interfaceInterpolator.solidInterfaceTemperature += self.solidTemperatureResidual

class AlgortihmBGSAitkenRelax(AlgortihmBGSStaticRelax):
    
    def __init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, Criterion, nbFSIIterMax, deltaT, totTime, timeIterTreshold=-1, omegaMax=1.0, mpiComm=None):
        """
        Des.
        """
        
        AlgortihmBGSStaticRelax.__init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, Criterion, nbFSIIterMax, deltaT, totTime, timeIterTreshold, omegaMax, mpiComm)
        
        ns = self.interfaceInterpolator.getNs()
        self.solidInterfaceResidualkM1 = FlexInterfaceData(ns, 3, self.mpiComm)
        self.aitkenCrit = 'max'
    
    def setOmega(self):
        """
        Des.
        """
        
        if self.FSIIter != 0:
            # --- Compute the dynamic Aitken coefficient --- #
            deltaInterfaceResidual = self.solidInterfaceResidual - self.solidInterfaceResidualkM1
            
            prodScalRes_X, prodScalRes_Y, prodScalRes_Z = deltaInterfaceResidual.dot(self.solidInterfaceResidualkM1)
            prodScalRes = prodScalRes_X + prodScalRes_Y + prodScalRes_Z
            
            deltaInterfaceResidual_NormX, deltaInterfaceResidual_NormY, deltaInterfaceResidual_NormZ = deltaInterfaceResidual.norm()
            deltaResNormSquare = deltaInterfaceResidual_NormX**2 + deltaInterfaceResidual_NormY**2 + deltaInterfaceResidual_NormZ**2
            
            if deltaResNormSquare != 0.:
                self.omega *= -prodScalRes/deltaResNormSquare
            else:
                self.omega = 0.1
        
        else:
            if self.aitkenCrit == 'max':
                self.omega = max(self.omegaMax, self.omega)
            else:
                self.omega = min(self.omegaMax, self.omega)
        
        self.omega = min(self.omega, 1.0)
        self.omega = max(self.omega, 0.1)
        
        mpiPrint('Aitken under-relaxation step with parameter {}'.format(self.omega), self.mpiComm)
        
        # --- Update the value of the residual for the next FSI iteration --- #
        self.solidInterfaceResidualkM1 = self.solidInterfaceResidual.copy()

class AlgortihmIQN_ILS(AlgortihmBGSAitkenRelax):
    """
    Des.
    """
    
    def __init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, Criterion, nbFSIIterMax, deltaT, totTime, timeIterTreshold=-1, omegaMax=1.0, nbTimeToKeep=0, computeTangentMatrixBasedOnFirstIt = False, mpiComm=None):
        """
        Des.
        """
        
        AlgortihmBGSAitkenRelax.__init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, Criterion, nbFSIIterMax, deltaT, totTime, timeIterTreshold, omegaMax, mpiComm)
        
        # --- Number of previous time steps used in the approximation of the tangent matrix --- #
        self.nbTimeToKeep = nbTimeToKeep
        
        # --- Option which allows to build the tangent matrix of a given time step using differences with respect to the first FSI iteration (delta_r_k = r_k+1 - r_0) instead of the previous iteration (delta_r_k = r_k+1 - r_k) --- #
        self.computeTangentMatrixBasedOnFirstIt = computeTangentMatrixBasedOnFirstIt
        
        # --- Option which determines the way the c coefficients are computes either using Degroote's QR decompoistion or simply using np.linalg.lstsq
        self.useQR = False
        self.tollQR = 1.0e-6 # tolerance used to get the tolerance for backward substitution after QR decomposition, toll, as toll = self.tollQR*norm(R)
        
        # --- Global V and W matrices for IQN-ILS algorithm, including information from previous time steps --- #
        self.V = []
        self.W = []

    def fsiCoupling(self, timeIter, time):
        """
        Interface Quasi Newton - Inverse Least Square (IQN-ILS) method for strong coupling FSI
        """
        
        if timeIter > self.timeIterTreshold:
            nbFSIIter = self.nbFSIIterMax
            mpiPrint('\n*************** Enter Interface Quasi Newton - Inverse Least Square (IQN-ILS) method for strong coupling FSI ***************', self.mpiComm)
        else:
             nbFSIIter = 1
        
        self.FSIIter = 0
        self.FSIConv = False
        self.errValue = 1.0
        self.errValue_CHT = 1.0 # Just for compatibility. CHT not yet implemented for the IQN-ILS algorithm.
        
        ns = self.interfaceInterpolator.getNs()
        
        # --- Initialize all the quantities used in the IQN-ILS method --- #
        res = FlexInterfaceData(ns, 3, self.mpiComm)
        solidInterfaceResidual0 = FlexInterfaceData(ns, 3, self.mpiComm)
        
        solidInterfaceDisplacement_tilde = FlexInterfaceData(ns, 3, self.mpiComm)
        solidInterfaceDisplacement_tilde1 = FlexInterfaceData(ns, 3, self.mpiComm)

        delta_ds = FlexInterfaceData(ns, 3, self.mpiComm)
        
        Vk_mat = np.zeros((3*ns,1))
        Wk_mat = np.zeros((3*ns,1))

        delta_ds_loc_X = np.zeros(0)
        delta_ds_loc_Y = np.zeros(0)
        delta_ds_loc_Z = np.zeros(0)
        
        if self.nbTimeToKeep!=0 and timeIter > self.nbTimeToKeep and not len(self.V) < self.nbTimeToKeep: # If information from previous time steps is re-used then Vk = V, Wk = W
            Vk = copy.deepcopy(self.V)
            Wk = copy.deepcopy(self.W)
            
        else: # If information from previous time steps is not re-used then Vk and Wk are empty lists of np.array()
            Vk = []
            Wk = []
        
        nIt = 0
        
        while ((self.FSIIter < nbFSIIter) and (not self.criterion.isVerified(self.errValue,self.errValue_CHT))):
            mpiPrint("\n>>>> FSI iteration {} <<<<\n".format(self.FSIIter), self.mpiComm)
            
            # --- Mesh morphing step (displacements interpolation, displacements communication, and mesh morpher call) --- #
            self.interfaceInterpolator.interpolateSolidDisplacementOnFluidMesh()
            self.interfaceInterpolator.setDisplacementToFluidSolver(time)
            mpiPrint('\nPerforming mesh deformation...\n', self.mpiComm)
            self.FluidSolver.meshUpdate(timeIter)
            
            # --- Fluid solver call for FSI subiteration --- #
            mpiPrint('\nLaunching fluid solver...', self.mpiComm)
            self.FluidSolver.run(time-self.deltaT, time)
            mpiBarrier(self.mpiComm)
            
            # --- Surface fluid loads interpolation and communication --- #
            mpiPrint('\nProcessing interface fluid loads...\n', self.mpiComm)
            self.interfaceInterpolator.getLoadsFromFluidSolver()
            mpiBarrier(self.mpiComm)
            if timeIter > self.timeIterTreshold:
                self.interfaceInterpolator.interpolateFluidLoadsOnSolidMesh()
                self.interfaceInterpolator.setLoadsToSolidSolver(time)
                
                # --- Solid solver call for FSI subiteration --- #
                mpiPrint('\nLaunching solid solver...\n', self.mpiComm)
                if self.myid in self.manager.getSolidSolverProcessors():
                    self.SolidSolver.run(time-self.deltaT, time)
                
                # --- Compute and monitor the FSI residual --- #
                res = self.computeSolidInterfaceResidual()
                self.errValue = self.criterion.update(res)
                mpiPrint('\nFSI error value : {}\n'.format(self.errValue), self.mpiComm)
                self.FSIConv = self.criterion.isVerified(self.errValue)
                
                # --- Initialize d_tilde for the construction of the Wk matrix -- #
                if self.myid in self.manager.getSolidInterfaceProcessors():
                    localSolidInterfaceDisp_X, localSolidInterfaceDisp_Y, localSolidInterfaceDisp_Z = self.SolidSolver.getNodalDisplacements()
                    for iVertex in range(self.manager.getNumberOfLocalSolidInterfaceNodes()):
                        iGlobalVertex = self.manager.getGlobalIndex('solid', self.myid, iVertex)
                        solidInterfaceDisplacement_tilde[iGlobalVertex] = [localSolidInterfaceDisp_X[iVertex], localSolidInterfaceDisp_Y[iVertex], localSolidInterfaceDisp_Z[iVertex]]
                
                solidInterfaceDisplacement_tilde.assemble()
                
                if self.FSIIter == 0 and (self.nbTimeToKeep == 0 or timeIter <= self.nbTimeToKeep or len(self.V) < self.nbTimeToKeep): # If information from previous time steps is re-used then this step is only performed at the first iteration of the first time step, otherwise it is performed at the first iteration of every time step
                    # --- Relax the solid position --- #
                    mpiPrint('\nProcessing interface displacements...\n', self.mpiComm)
                    self.relaxSolidPosition()
                
                else:
                    # --- Construct Vk and Wk matrices for the computation of the approximated tangent matrix --- #
                    mpiPrint('\nCorrect solid interface displacements using IQN-ILS method...\n', self.mpiComm)
                    
                    # --- Start gathering on root process --- #
                    res_X_Gat, res_Y_Gat, res_Z_Gat = mpiGatherInterfaceData(res, ns, self.mpiComm, 0)
                    solidInterfaceResidual0_X_Gat, solidInterfaceResidual0_Y_Gat, solidInterfaceResidual0_Z_Gat = mpiGatherInterfaceData(solidInterfaceResidual0, ns, self.mpiComm, 0)
                    solidInterfaceDisplacement_tilde_X_Gat, solidInterfaceDisplacement_tilde_Y_Gat, solidInterfaceDisplacement_tilde_Z_Gat = mpiGatherInterfaceData(solidInterfaceDisplacement_tilde, ns, self.mpiComm, 0)
                    solidInterfaceDisplacement_tilde1_X_Gat, solidInterfaceDisplacement_tilde1_Y_Gat, solidInterfaceDisplacement_tilde1_Z_Gat = mpiGatherInterfaceData(solidInterfaceDisplacement_tilde1, ns, self.mpiComm, 0)
                    
                    if self.myid == 0:
                        if self.FSIIter > 0: # Either information from previous time steps is re-used or not, Vk and Wk matrices are enriched only starting from the second iteration of every FSI loop
                            delta_res = np.concatenate([res_X_Gat - solidInterfaceResidual0_X_Gat, res_Y_Gat - solidInterfaceResidual0_Y_Gat, res_Z_Gat - solidInterfaceResidual0_Z_Gat], axis=0)
                            delta_d = np.concatenate([solidInterfaceDisplacement_tilde_X_Gat - solidInterfaceDisplacement_tilde1_X_Gat, solidInterfaceDisplacement_tilde_Y_Gat - solidInterfaceDisplacement_tilde1_Y_Gat, solidInterfaceDisplacement_tilde_Z_Gat - solidInterfaceDisplacement_tilde1_Z_Gat], axis = 0)
                        
                            Vk.insert(0, delta_res)
                            Wk.insert(0, delta_d)
                        
                            nIt+=1
                    
                        Vk_mat = np.vstack(Vk).T
                        Wk_mat = np.vstack(Wk).T
                        
                        if self.useQR: # Technique described by Degroote et al.
                            Q, R = sp.linalg.qr(Vk_mat, mode='economic')
                            
                            s = np.dot(np.transpose(Q), -np.concatenate([res_X_Gat, res_Y_Gat, res_Z_Gat], axis=0))
                            
                            toll = self.tollQR*sp.linalg.norm(R, 2)
                            c = solve_upper_triangular_mod(R, s, toll)
                            
                            delta_ds_loc = np.split((np.dot(Wk_mat,c).T + np.concatenate([res_X_Gat, res_Y_Gat, res_Z_Gat], axis=0)),3,axis=1)

                            delta_ds_loc_X = np.concatenate(delta_ds_loc[0])
                            delta_ds_loc_Y = np.concatenate(delta_ds_loc[1])
                            delta_ds_loc_Z = np.concatenate(delta_ds_loc[2])
                        else:
                            c = np.linalg.lstsq(Vk_mat, -np.concatenate([res_X_Gat, res_Y_Gat, res_Z_Gat], axis=0))[0]
                            
                            delta_ds_loc = np.split((np.dot(Wk_mat,c).T + np.concatenate([res_X_Gat, res_Y_Gat, res_Z_Gat], axis=0)),3,axis=0)
                            
                            delta_ds_loc_X = delta_ds_loc[0]
                            delta_ds_loc_Y = delta_ds_loc[1]
                            delta_ds_loc_Z = delta_ds_loc[2]
                        
                        for iVertex in range(delta_ds_loc_X.shape[0]):
                            iGlobalVertex = self.manager.getGlobalIndex('solid', self.myid, iVertex)
                            delta_ds[iGlobalVertex] = [delta_ds_loc_X[iVertex], delta_ds_loc_Y[iVertex], delta_ds_loc_Z[iVertex]]
                    
                    # --- Go back to parallel run --- #
                    mpiBarrier(self.mpiComm)
                    delta_ds.assemble()
                    self.interfaceInterpolator.solidInterfaceDisplacement += delta_ds
                    
                if self.computeTangentMatrixBasedOnFirstIt:
                    if self.FSIIter == 0:
                        solidInterfaceResidual0 = res.copy()
                        solidInterfaceDisplacement_tilde1 = solidInterfaceDisplacement_tilde.copy()
                else:
                    solidInterfaceResidual0 = res.copy()
                    solidInterfaceDisplacement_tilde1 = solidInterfaceDisplacement_tilde.copy()
            
            if self.writeInFSIloop == True:
                self.writeRealTimeData(timeIter, time, self.errValue, 0.0)
            
            self.FSIIter += 1
        
        # if comm.myself == rootProcess
        # update of the matrices V and W at the end of the while
        if self.nbTimeToKeep != 0 and timeIter > self.timeIterTreshold and self.FSIIter > 1:
            
            if (self.FSIIter >= nbFSIIter):
                mpiPrint('WARNING: IQN-ILS using information from {} previous time steps reached max number of iterations. Next time step is run without using any information from previous time steps!'.format(self.nbTimeToKeep), self.mpiComm)
                
                self.V = []
                self.W = []
                
            else:
                mpiPrint('\nUpdating V and W matrices...\n', self.mpiComm)
                
                self.V.insert(0, Vk_mat[:,0:nIt].T)
                self.W.insert(0, Wk_mat[:,0:nIt].T)
                
                if timeIter > self.nbTimeToKeep and not len(self.V) <= self.nbTimeToKeep:
                    del self.V[-1]
                    del self.W[-1]
            
        # --- Update the FSI history file --- #
        if timeIter > self.timeIterTreshold:
            mpiPrint('\n*************** IQN-ILS is converged ***************', self.mpiComm)

# --- Solid test algorithm ---
class FsiSolidTestAlgorithm:
    def __init__(self, _solid):
        self.solid = _solid
        
    def run(self):
        # --------------------------
        # fake FSI solver
        # --------------------------
        
        t1 = 0.0  # initial time
        dt = 0.5  # time step size
        nt = 10
        
        # we want nt time steps
        for j in range(nt):
            
            # each time step is arbitrarily calculated twice (for testing purpose)
            for i in range(2):
                
                t2=t1+dt  # time to be calculated
                
                self.solid.fakeFluidSolver(t2)  # creates some dummy loads for time=t2
                
                # run solid solver
                print '='*80
                print "running from %f to %f: try #%d" % (t1,t2,i+1)
                print '='*80
                self.solid.run(t1,t2)
                
                # gets the deformed interface
                dx, dy, dz = self.solid.getNodalDisplacements()
                print dx
                print dy
                print dz
                
            self.solid.update()
            self.solid.save()
            
            t1=t2 # fsi loop has converged - time t2 is accepted
        
        # end.
