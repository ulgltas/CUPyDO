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
import os, os.path, sys, time, string

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
        #interfData_X_Gat = mpiGatherv(interfData.getXArray(), localSize, globalSize, mpiComm, 0)
        #interfData_Y_Gat = mpiGatherv(interfData.getYArray(), localSize, globalSize, mpiComm, 0)
        #interfData_Z_Gat = mpiGatherv(interfData.getZArray(), localSize, globalSize, mpiComm, 0)
        #return (interfData_X_Gat, interfData_Y_Gat, interfData_Z_Gat)
    else:
        for iDim in range(interfData.nDim):
            interfData_Gat.append(interfData.getData(iDim))
        #return (interfData.getDataX(), interfData.getDataY(), interfData.getDataZ())

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
                newData.dataContainer[iDim].copy(self.dataContainer[iDim])

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
        
        # --- Create the array for external communication (displacement, velocity and velocity at the previous time step --- #
        self.nodalDisp_X = np.zeros(self.nPhysicalNodes)
        self.nodalDisp_Y = np.zeros(self.nPhysicalNodes)
        self.nodalDisp_Z = np.zeros(self.nPhysicalNodes)
        self.nodalVel_X = np.zeros(self.nPhysicalNodes)
        self.nodalVel_Y = np.zeros(self.nPhysicalNodes)
        self.nodalVel_Z = np.zeros(self.nPhysicalNodes)
        self.nodalVel_XNm1 = np.zeros(self.nPhysicalNodes)
        self.nodalVel_YNm1 = np.zeros(self.nPhysicalNodes)
        self.nodalVel_ZNm1 = np.zeros(self.nPhysicalNodes)
    
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
    
    def setInitialMeshDeformation(self):
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
    
    def applyNodalDisplacements(self, dx, dy, dz, dx_nM1, dy_nM1, dz_nM1, haloNodesDisplacements,time):
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
    
    def __init__(self, tolerance):
        """
        Description.
        """
        
        self.tol = tolerance
        self.epsilon = 0
    
    def isVerified(self, epsilon):
        """
        Description.
        """
        
        if epsilon < self.tol:
            return True
        else:
            return False

class DispNormCriterion(Criterion):
    """
    Description.
    """
    
    def __init__(self, tolerance):
        """
        Description.
        """
        
        self.tol = tolerance
    
    def update(self, res):
        """
        Des.
        """
        
        normX, normY, normZ = res.norm()
        
        norm = sqrt(normX**2 + normY**2 + normZ**2)
        
        self.epsilon = norm
        
        return self.epsilon

# ----------------------------------------------------------------------
#    Manager class
# ----------------------------------------------------------------------

class Manager:
    """
    Description.
    """

    def __init__(self, FluidSolver, SolidSolver, nDim, computationType='steady', withCht=False, mpiComm=None):
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
        self.withCht = withCht

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

    def getFluidIndexing(self):
        """
        Des.
        """

        return self.fluidIndexing

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

    def __init__(self, Manager, FluidSolver, SolidSolver, chtTransferMethod = 'TFFB', mpiComm = None):
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
        if self.manager.withCht :
            self.chtTransferMethod = chtTransferMethod
        #FFTB = Flux Forward Temperature Back
        #TFFB = Temperature Forward Flux Back
        #hFTB = Heat Transfer Coefficient Forward Temperature Back
        #hFFB = Heat Transfer Coefficient Forward Flux Back

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

        #self.solidInterfaceHeatFlux = None
        #self.fluidInterfaceHeatFlux = None
        #self.solidInterfaceTemperature = None
        #self.fluidInterfaceTemperature = None
        #if self.managet.withCht :
        #    self.solidInterfaceHeatFlux = InterfaceData(self.ns + self.d, self.mpiComm)
        #    self.fluidInterfaceHeatFlux = InterfaceData(self.nf, self.mpiComm)
        #    self.solidInterfaceTemperature = InterfaceData(self.ns + self.d, self.mpiComm)
        #    self.fluidInterfaceTemperature = InterfaceData(self.nf, self.mpiComm)

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

    def redistributeDataToFluidSolver(self, fluidInterfaceData):
        """
        Description
        """

        localFluidInterfaceData_array_X = None
        localFluidInterfaceData_array_Y = None
        localFluidInterfaceData_array_Z = None

        if self.mpiComm != None:
            localSize = fluidInterfaceData.getXArray().shape[0]
            fluidInterfaceData_array_X_recon = mpiGatherv(self.fluidInterfaceData.getXArray(), localSize, self.nf, self.mpiComm, 0)
            fluidInterfaceData_array_Y_recon = mpiGatherv(self.fluidInterfaceData.getYArray(), localSize, self.nf, self.mpiComm, 0)
            fluidInterfaceData_array_Z_recon = mpiGatherv(self.fluidInterfaceData.getZArray(), localSize, self.nf, self.mpiComm, 0)
            haloNodesData = {}
            if self.myid == 0:
                for iProc in self.manager.getFluidInterfaceProcessors():
                    fluidPhysicalInterfaceNodesDistribution = self.manager.getFluidPhysicalInterfaceNodesDistribution()
                    fluidGlobalIndexRange = self.manager.getFluidGlobalIndexRange()
                    sendBuff_X = np.zeros(fluidPhysicalInterfaceNodesDistribution[iProc])
                    sendBuff_Y = np.zeros(fluidPhysicalInterfaceNodesDistribution[iProc])
                    sendBuff_Z = np.zeros(fluidPhysicalInterfaceNodesDistribution[iProc])
                    globalIndex = fluidGlobalIndexRange[iProc][iProc][0]
                    sendBuffHalo = {}
                    for iVertex in range(fluidPhysicalInterfaceNodesDistribution[iProc]):
                        sendBuff_X[iVertex] = fluidInterfaceData_array_X_recon[globalIndex]
                        sendBuff_Y[iVertex] = fluidInterfaceData_array_Y_recon[globalIndex]
                        sendBuff_Z[iVertex] = fluidInterfaceData_array_Z_recon[globalIndex]
                        globalIndex += 1
                    fluidHaloNodesList = self.manager.getFluidHaloNodesList()
                    fluidIndexing = self.manager.getFluidIndexing()
                    for key in fluidHaloNodesList[iProc].keys():
                        globalIndex = fluidIndexing[key]
                        sendBuffHalo[key] = (fluidInterfaceData_array_X_recon[globalIndex], fluidInterfaceData_array_Y_recon[globalIndex], fluidInterfaceData_array_Z_recon[globalIndex])
                    self.mpiComm.Send(sendBuff_X, dest=iProc, tag = 1)
                    self.mpiComm.Send(sendBuff_Y, dest=iProc, tag = 2)
                    self.mpiComm.Send(sendBuff_Z, dest=iProc, tag = 3)
                    self.mpiComm.send(sendBuffHalo, dest = iProc, tag=4)
            if self.myid in self.manager.getFluidInterfaceProcessors():
                localFluidInterfaceData_array_X = np.zeros(self.nf_loc)
                localFluidInterfaceData_array_Y = np.zeros(self.nf_loc)
                localFluidInterfaceData_array_Z = np.zeros(self.nf_loc)
                self.mpiComm.Recv(localFluidInterfaceData_array_X, source=0, tag=1)
                self.mpiComm.Recv(localFluidInterfaceData_array_Y, source=0, tag=2)
                self.mpiComm.Recv(localFluidInterfaceData_array_Z, source=0, tag=3)
                haloNodesData = self.mpiComm.recv(source=0, tag=4)

        return (localFluidInterfaceData_array_X, localFluidInterfaceData_array_Y, localFluidInterfaceData_array_Z, haloNodesData)

    def setLoadsToSolidSolver(self, time):
        """
        Des.
        """

        #self.checkTotalLoad()

        FFX, FFY, FFZ = self.fluidInterfaceLoads.sum()
        FX = 0.0
        FY = 0.0
        FZ = 0.0
        
        FXT = 0.
        FYT = 0.
        FZT = 0.
        
        # --- Redistribute the interpolated solid loads according to the partitions that own the solid interface --- #

        if self.mpiComm != None:
            localSize = self.solidInterfaceLoads.getDataArray(0).shape[0]
            solidLoads_array_recon = []
            for iDim in range(self.solidInterfaceLoads.nDim):
                array_recon = mpiGatherv(self.solidInterfaceLoads.getDataArray(iDim), localSize, self.ns+self.d, self.mpiComm, 0)
                solidLoads_array_recon.append(array_recon)

            if self.myid == 0:
                for iProc in self.manager.getSolidInterfaceProcessors():
                    solidPhysicalInterfaceNodesDistribution = self.manager.getSolidPhysicalInterfaceNodesDistribution()
                    solidGlobalIndexRange = self.manager.getSolidGlobalIndexRange()
                    sendBuff = []
                    for iDim in range(self.solidInterfaceLoads.nDim):
                        sendBuff_i = np.zeros(solidPhysicalInterfaceNodesDistribution[iProc])
                        sendBuff.append(sendBuff_i)
                    globalIndex = solidGlobalIndexRange[iProc][iProc][0]
                    for iVertex in range(solidPhysicalInterfaceNodesDistribution[iProc]):
                        for iDim in range(self.solidInterfaceLoads.nDim):
                            sendBuff[iDim][iVertex] = solidLoads_array_recon[iDim][globalIndex]
                        globalIndex += 1
                    iTagSend = 1
                    for iDim in range(self.solidInterfaceLoads.nDim):
                        self.mpiComm.Send(sendBuff[iDim], dest=iProc, tag = iTagSend)
                        iTagSend += 1
            if self.myid in self.manager.getSolidInterfaceProcessors():
                localSolidLoads_array = []
                iTagRec = 1
                for iDim in range(self.solidInterfaceLoads.nDim):
                    local_array = np.zeros(self.ns_loc)
                    self.mpiComm.Recv(local_array, source=0, tag = iTagRec)
                    localSolidLoads_array.append(local_array)
                    iTagRec += 1
                self.SolidSolver.applyNodalLoads(localSolidLoads_array[0], localSolidLoads_array[1], localSolidLoads_array[2], time)
                for iVertex in range(self.ns_loc):
                    FX += localSolidLoads_array[0][iVertex]
                    FY += localSolidLoads_array[1][iVertex]
                    FZ += localSolidLoads_array[2][iVertex]
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
            localSize = self.fluidInterfaceDisplacement.getDataArray(0).shape[0]
            fluidInterface_array_Disp_recon = []
            for iDim in range(self.fluidInterfaceDisplacement.nDim):
                array_recon = mpiGatherv(self.fluidInterfaceDisplacement.getDataArray(iDim), localSize, self.nf, self.mpiComm, 0)
                fluidInterface_array_Disp_recon.append(array_recon)
            haloNodesDisplacements = {}
            if self.myid == 0:
                for iProc in self.manager.getFluidInterfaceProcessors():
                    fluidPhysicalInterfaceNodesDistribution = self.manager.getFluidPhysicalInterfaceNodesDistribution()
                    fluidGlobalIndexRange = self.manager.getFluidGlobalIndexRange()
                    sendBuff = []
                    for iDim in range(self.fluidInterfaceDisplacement.nDim):
                        sendBuff_i = np.zeros(fluidPhysicalInterfaceNodesDistribution[iProc])
                        sendBuff.append(sendBuff_i)
                    globalIndex = fluidGlobalIndexRange[iProc][iProc][0]
                    sendBuffHalo = {}
                    for iVertex in range(fluidPhysicalInterfaceNodesDistribution[iProc]):
                        for iDim in range(self.fluidInterfaceDisplacement.nDim):
                            sendBuff[iDim][iVertex] = fluidInterface_array_Disp_recon[iDim][globalIndex]
                        globalIndex += 1
                    fluidHaloNodesList = self.manager.getFluidHaloNodesList()
                    fluidIndexing = self.manager.getFluidIndexing()
                    for key in fluidHaloNodesList[iProc].keys():
                        globalIndex = fluidIndexing[key]
                        sendBuffHalo[key] = (fluidInterface_array_Disp_recon[0][globalIndex], fluidInterface_array_Disp_recon[1][globalIndex], fluidInterface_array_Disp_recon[2][globalIndex])
                    iTagSend = 1
                    for iDim in range(self.fluidInterfaceDisplacement.nDim):
                        self.mpiComm.Send(sendBuff[iDim], dest=iProc, tag = iTagSend)
                        iTagSend += 1
                    self.mpiComm.send(sendBuffHalo, dest = iProc, tag=iTagSend)
            if self.myid in self.manager.getFluidInterfaceProcessors():
                localFluidInterface_array_Disp = []
                iTagRec = 1
                for iDim in range(self.fluidInterfaceDisplacement.nDim):
                    local_array = np.zeros(self.nf_loc)
                    self.mpiComm.Recv(local_array, source=0, tag=iTagRec)
                    localFluidInterface_array_Disp.append(local_array)
                    iTagRec += 1
                haloNodesDisplacements = self.mpiComm.recv(source=0, tag=iTagRec)
                self.FluidSolver.applyNodalDisplacements(localFluidInterface_array_Disp[0], localFluidInterface_array_Disp[1], localFluidInterface_array_Disp[2], localFluidInterface_array_Disp[0], localFluidInterface_array_Disp[1], localFluidInterface_array_Disp[2], haloNodesDisplacements, time)
        else:
            self.FluidSolver.applyNodalDisplacements(self.fluidInterfaceDisplacement.getDataArray(0), self.fluidInterfaceDisplacement.getDataArray(1), self.fluidInterfaceDisplacement.getDataArray(2), self.fluidInterfaceDisplacement.getDataArray(0), self.fluidInterfaceDisplacement.getDataArray(1), self.fluidInterfaceDisplacement.getDataArray(2), {}, time)

    def setHeatFluxToFluidSolver(self, time):
        """
        Description.
        """

        if self.mpiComm != None:
            (localFluidInterface_array_QX, localFluidInterface_array_QY, localFluidInterface_array_QZ, haloNodesHeatFlux) = self.redistributeDataToFluidSolver(self.fluidInterfaceHeatFlux)
            if self.myid in self.manager.getFluidInterfaceProcessors():
                self.FluidSolver.applyNodalHeatFluxes(localFluidInterface_array_QX, localFluidInterface_array_QY, localFluidInterface_array_QZ, localFluidInterface_array_QX, localFluidInterface_array_QY, localFluidInterface_array_QZ, haloNodesHeatFlux, time)
        else:
            self.FluidSolver.applyNodalHeatFluxes(self.fluidInterfaceHeatFlux.getXArray(), self.fluidInterfaceHeatFlux.getYArray(), self.fluidInterfaceHeatFlux.getZArray(), self.fluidInterfaceHeatFlux.getXArray(), self.fluidInterfaceHeatFlux.getYArray(), self.fluidInterfaceHeatFlux.getZArray(), {}, time)

    def setTemperatureToSolidSolver(self, time):
        """
        Description
        """

        if self.mpiComm != None:
            (localSolidInterface_array_TX, localSolidInterface_array_QY, localSolidInterface_array_QZ, haloNodesHeatFlux)

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

    def __init__(self, Manager, FluidSolver, SolidSolver, chtTransferMethod = 'TFFB', mpiComm = None):
        """
        Description
        """

        InterfaceInterpolator.__init__(self, Manager, FluidSolver, SolidSolver, chtTransferMethod, mpiComm)

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

        self.H.mult(self.solidInterfaceHeatFlux.getDataX(), self.fluidInterfaceHeatFlux.getDataX())
        self.H.mult(self.solidInterfaceHeatFlux.getDataY(), self.fluidInterfaceHeatFlux.getDataY())
        self.H.mult(self.solidInterfaceHeatFlux.getDataZ(), self.fluidInterfaceHeatFlux.getDataZ())

    def interpolateSolidTemperatureOnFluidMesh(self):
        """
        Description
        """

        self.H.mult(self.solidInterfaceTemperature.getDataX(), self.fluidInterfaceTemperature.getDataX())
        self.H.mult(self.solidInterfaceTemperature.getDataY(), self.fluidInterfaceTemperature.getDataY())
        self.H.mult(self.solidInterfaceTemperature.getDataZ(), self.fluidInterfaceTemperature.getDataZ())

    def interpolateFluidHeatFluxOnSolidMesh(self):
        """
        Description.
        """

        self.H_T.mult(self.fluidInterfaceHeatFlux.getDataX(), self.solidInterfaceHeatFlux.getDataX())
        self.H_T.mult(self.fluidInterfaceHeatFlux.getDataY(), self.solidInterfaceHeatFlux.getDataY())
        self.H_T.mult(self.fluidInterfaceHeatFlux.getDataZ(), self.solidInterfaceHeatFlux.getDataZ())

    def interpolateFluidTemperatureOnSolidMesh(self):
        """
        Description.
        """

        self.H_T.mult(self.fluidInterfaceTemperature.getDataX(), self.solidInterfaceTemperature.getDataX())
        self.H_T.mult(self.fluidInterfaceTemperature.getDataY(), self.solidInterfaceTemperature.getDataY())
        self.H_T.mult(self.fluidInterfaceTemperature.getDataZ(), self.solidInterfaceTemperature.getDataZ())
    

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
        else:
            self.interfaceInterpolator.getDisplacementFromSolidSolver()
    
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
        self.solidInterfaceResidual = predictedDisplacement - self.interfaceInterpolator.solidInterfaceDisplacement
        
        return self.solidInterfaceResidual
    
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
        histFile.write("TimeIter\tTime\tFSIError\tFSINbIter\n")
        histFile.close()
    
    def writeRealTimeData(self, timeIter, time, error):
        """
        Des
        """
        
        if self.myid == 0:
            self.FluidSolver.saveRealTimeData(time, self.FSIIter)
            if timeIter >= self.timeIterTreshold:
                self.SolidSolver.saveRealTimeData(time, self.FSIIter)
            histFile = open('FSIhistory.ascii', "a")
            histFile.write(str(timeIter) + '\t' + str(time) + '\t' + str(error) + '\t' + str(self.FSIIter) + '\n')
            histFile.close()
    
    def run(self):
        """
        Des.
        """
        
        # --- Initialize output manager --- #
        self.iniRealTimeData()
        
        mpiPrint('\n**********************************', self.mpiComm)
        mpiPrint('*         Begin FSI computation            *', self.mpiComm)
        mpiPrint('**********************************\n', self.mpiComm)
        
        #If no restart
        mpiPrint('Setting FSI initial conditions...', self.mpiComm)
        self.setFSIInitialConditions(0.0)
        mpiPrint('\nFSI initial conditions are set', self.mpiComm)
        
        if self.manager.computationType == 'unsteady':
            self.__unsteadyRun()
        else:
            time = self.totTime
            timeIter = 0
            self.deltaT = self.totTime
            self.writeInFSIloop = True
            self.fsiCoupling(timeIter, time)
        
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
            
            # --- Update the fluid and solid solver for the next time step --- #
            if self.myid in self.manager.getSolidSolverProcessors():
                self.SolidSolver.update()
            self.FluidSolver.update(self.deltaT)
            
            # --- Write fluid and solid solution, update FSI history  ---#
            self.FluidSolver.save(timeIter)
            
            if self.myid in self.manager.getSolidSolverProcessors():
                self.SolidSolver.save()
            
            self.writeRealTimeData(timeIter, time, self.errValue)
            
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
        
        while ((self.FSIIter < nbFSIIter) and (not self.criterion.isVerified(self.errValue))):
            mpiPrint("\n>>>> FSI iteration {} <<<<\n".format(self.FSIIter), self.mpiComm)
            
            # --- Mesh morphing step (displacements interpolation, displacements communication, and mesh morpher call) --- #
            self.interfaceInterpolator.interpolateSolidDisplacementOnFluidMesh()
            self.interfaceInterpolator.setDisplacementToFluidSolver(time)
            mpiPrint('\nPerforming mesh deformation...\n', self.mpiComm)
            self.FluidSolver.meshUpdate(timeIter)
            #self.interpolator.chtTransferMethod = 'TFFB':
                #self.interpolateSolidHeatFluxOnFluidMesh()
                #self.setHeatFluxToFluidSolver(time)
            
            # --- Fluid solver call for FSI subiteration --- #
            mpiPrint('\nLaunching fluid solver...', self.mpiComm)
            self.FluidSolver.run(time-self.deltaT, time)
            mpiBarrier(self.mpiComm)
            
            # --- Surface fluid loads interpolation and communication --- #
            mpiPrint('\nProcessing interface fluid loads...\n', self.mpiComm)
            self.interfaceInterpolator.getLoadsFromFluidSolver()
            #self.interpolator.chtTransferMethod = 'TFFB':
                #self.interpolator.getTemperatureFromFluidSolver()
            mpiBarrier(self.mpiComm)
            if timeIter > self.timeIterTreshold:
                self.interfaceInterpolator.interpolateFluidLoadsOnSolidMesh()
                self.interfaceInterpolator.setLoadsToSolidSolver(time)
                #self.interpolator.chtTransferMethod = 'TFFB':
                  #self.interpolator.interpolateFluidTemperatureOnSolidMesh
                  #self.interpolator.setTemperatureToSolidSolver(time)
                
                # --- Solid solver call for FSI subiteration --- #
                mpiPrint('\nLaunching solid solver...\n', self.mpiComm)
                if self.myid in self.manager.getSolidSolverProcessors():
                    self.SolidSolver.run(time-self.deltaT, time)

                # --- Compute and monitor the FSI residual --- #
                res = self.computeSolidInterfaceResidual()
                self.errValue = self.criterion.update(res)
                mpiPrint('\nFSI error value : {}\n'.format(self.errValue), self.mpiComm)
                self.FSIConv = self.criterion.isVerified(self.errValue)
                
                # --- Relaxe the solid position --- #
                mpiPrint('\nProcessing interface displacements...\n', self.mpiComm)
                self.relaxSolidPosition()
            
            if self.writeInFSIloop == True:
                self.writeRealTimeData(timeIter, time, self.errValue)
            
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

class AlgortihmBGSAitkenRelax(AlgortihmBGSStaticRelax):
    
    def __init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, Criterion, nbFSIIterMax, deltaT, totTime, timeIterTreshold=-1, omegaMax=1.0, mpiComm=None):
        """
        Des.
        """
        
        AlgortihmBGSStaticRelax.__init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, Criterion, nbFSIIterMax, deltaT, totTime, timeIterTreshold, omegaMax, mpiComm)
        
        ns = self.interfaceInterpolator.getNs()
        self.solidInterfaceResidualkM1 = FlexInterfaceData(ns, 3, self.mpiComm)
    
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
            
            self.omega *= -prodScalRes/deltaResNormSquare
        else:
            self.omega = max(self.omegaMax, self.omega)
        
        self.omega = min(self.omega, 1.0)
        self.omega = max(self.omega, 0.0)
        
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
        
        while ((self.FSIIter < nbFSIIter) and (not self.criterion.isVerified(self.errValue))):
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
                        
                        Q, R = sp.linalg.qr(Vk_mat, mode='economic')
                        
                        s = np.dot(np.transpose(Q), -np.concatenate([res_X_Gat, res_Y_Gat, res_Z_Gat], axis=0))
                    
                        toll = 1e-10*np.linalg.norm(np.concatenate([res_X_Gat, res_Y_Gat, res_Z_Gat]))
                        c = solve_upper_triangular_mod(R, s, toll)
                    
                        delta_ds_loc = np.split((np.dot(Wk_mat,c).T + np.concatenate([res_X_Gat, res_Y_Gat, res_Z_Gat], axis=0)),3,axis=1)

                        delta_ds_loc_X = np.concatenate(delta_ds_loc[0])
                        delta_ds_loc_Y = np.concatenate(delta_ds_loc[1])
                        delta_ds_loc_Z = np.concatenate(delta_ds_loc[2])

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
                self.writeRealTimeData(timeIter, time, self.errValue)
            
            self.FSIIter += 1
        
        # if comm.myself == rootProcess
        # update of the matrices V and W at the end of the while
        if self.nbTimeToKeep != 0 and timeIter > self.timeIterTreshold and self.FSIIter > 1:
            
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
