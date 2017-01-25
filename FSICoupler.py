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
import os, os.path, sys, time, string

import socket, fnmatch
import fsi_pyutils

# global vars (underscore prevent them to be imported with "from module import *")
_theModule  = None
_theWDir    = None # workspace directory
_theWDirRoot = os.getcwd()  # base directory du workspace

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


# ----------------------------------------------------------------------
#    InterfaceData class
# ----------------------------------------------------------------------

class InterfaceData:
    """
    Des.
    """

    def __init__(self, nPoint, mpiComm=None):
        """
        Des.
        """

        self.mpiComm = mpiComm
        self.nPoint = nPoint

        if mpiComm != None:
            from petsc4py import PETSc
            self.data_X = PETSc.Vec().create(self.mpiComm)
            self.data_Y = PETSc.Vec().create(self.mpiComm)
            self.data_Z = PETSc.Vec().create(self.mpiComm)
            self.data_X.setType('mpi')
            self.data_Y.setType('mpi')
            self.data_Z.setType('mpi')
            self.data_X.setSizes(nPoint)
            self.data_Y.setSizes(nPoint)
            self.data_Z.setSizes(nPoint)
            self.myid = self.mpiComm.Get_rank()
            self.mpiSize = self.mpiComm.Get_size()
            startIndex , stopIndex = self.data_X.getOwnershipRange()
            self.indexing = self.mpiComm.allgather((startIndex , stopIndex))
        else:
            self.data_X = np.zeros(nPoint, dtype=float)
            self.data_Y = np.zeros(nPoint, dtype=float)
            self.data_Z = np.zeros(nPoint, dtype=float)
            self.myid = 0
            self.mpiSize = 1

    def setXValues(self, indices_list, values_array):
        """
        Des.
        """
        
        if self.mpiComm != None:
            self.data_X.setValues(indices_list, values_array)
        else:
            self.data_X[indices_list] = values_array

    def setYValues(self, indices_list, values_array):
        """
        Des.
        """
        
        if self.mpiComm != None:
            self.data_Y.setValues(indices_list, values_array)
        else:
            self.data_Y[indices_list] = values_array

    def setZValues(self, indices_list, values_array):
        """
        Des.
        """
        
        if self.mpiComm != None:
            self.data_Z.setValues(indices_list, values_array)
        else:
            self.data_Z[indices_list] = values_array

    def getDataX(self):
        """
        des.
        """

        return self.data_X

    def getDataY(self):
        """
        des.
        """

        return self.data_Y

    def getDataZ(self):
        """
        des.
        """

        return self.data_Z

    def getXArray(self):
        """
        Des.
        """

        if self.mpiComm != None:
            return self.data_X.getArray()
        else:
            return self.data_X

    def getYArray(self):
        """
        Des.
        """

        if self.mpiComm != None:
            return self.data_Y.getArray()
        else:
            return self.data_Y

    def getZArray(self):
        """
        Des.
        """

        if self.mpiComm != None:
            return self.data_Z.getArray()
        else:
            return self.data_Z

    def assemble(self):
        """
        des.
        """

        if self.mpiComm != None:
            self.data_X.assemblyBegin()
            self.data_X.assemblyEnd()
            self.data_Y.assemblyBegin()
            self.data_Y.assemblyEnd()
            self.data_Z.assemblyBegin()
            self.data_Z.assemblyEnd()

    def view(self):
        """
        Des.
        """

        self.data_X.view()

    def norm(self):
        """
        Des
        """

        if self.mpiComm != None:
            norm_X = self.data_X.norm()
            norm_Y = self.data_Y.norm()
            norm_Z = self.data_Z.norm()
        else:
            norm_X = np.linalg.norm(self.data_X)
            norm_Y = np.linalg.norm(self.data_Y)
            norm_Z = np.linalg.norm(self.data_Z)

        return (norm_X, norm_Y, norm_Z)

    def sum(self):
        """
        des.
        """

        sum_X = self.data_X.sum()
        sum_Y = self.data_Y.sum()
        sum_Z = self.data_Z.sum()

        return (sum_X, sum_Y, sum_Z)

    def copy(self):
        """
        Des.
        """

        newData = InterfaceData(self.nPoint, self.mpiComm)
        if self.mpiComm != None:
            self.data_X.copy(newData.data_X)
            self.data_Y.copy(newData.data_Y)
            self.data_Z.copy(newData.data_Z)
        else:
            newData.data_X = self.data_X.copy()
            newData.data_Y = self.data_Y.copy()
            newData.data_Z = self.data_Z.copy()

        return newData

    def dot(self, data):
        """
        Des.
        """

        if type(data) != type(self):
            #raise type error
            print "WRONG TYPE"
        else:
            dotX = self.data_X.dot(data.data_X)
            dotY = self.data_Y.dot(data.data_Y)
            dotZ = self.data_Z.dot(data.data_Z)

        return (dotX, dotY, dotZ)

    def __setitem__(self, index, values):
        """
        des.
        """

        if type(values) != tuple:
            #raise TypeError
            print "NEED TUPLE"
        else:
            val_x, val_y, val_z = values
            self.data_X[index] = val_x
            self.data_Y[index] = val_y
            self.data_Z[index] = val_z

    def __getitem__(self, index):
        """
        des.
        """

        if self.mpiComm != None:
            send = None
            rcv = None
            for iProc in range(self.mpiSize):
                start, stop = self.indexing[iProc]
                if index in range(start, stop):
                    sender = iProc
                    break
            mpiBarrier(self.mpiComm)
            if self.myid == sender:
                send = (float(self.data_X[index]), float(self.data_Y[index]), float(self.data_Z[index]))
            rcv = self.mpiComm.bcast(send, sender)
            return rcv
        else:
            return (float(self.data_X[index]), float(self.data_Y[index]), float(self.data_Z[index]))

    def __add__(self, dataToAdd):
        """
        Des.
        """

        newData = InterfaceData(self.nPoint, self.mpiComm)
        if type(dataToAdd) == type(self):
            newData.data_X = self.data_X + dataToAdd.data_X
            newData.data_Y = self.data_Y + dataToAdd.data_Y
            newData.data_Z = self.data_Z + dataToAdd.data_Z
        elif type(dataToAdd) == float:
            newData.data_X = self.data_X + dataToAdd
            newData.data_Y = self.data_Y + dataToAdd
            newData.data_Z = self.data_Z + dataToAdd
        elif type(dataToAdd) == int:
            newData.data_X = self.data_X + float(dataToAdd)
            newData.data_Y = self.data_Y + float(dataToAdd)
            newData.data_Z = self.data_Z + float(dataToAdd)

        return newData

    def __radd__(self, dataToAdd):
        """
        des.
        """

        newData = self + dataToAdd

        return newData

    def __iadd__(self, dataToAdd):
        """
        Des.
        """

        if type(dataToAdd) == type(self):
            self.data_X += dataToAdd.data_X
            self.data_Y += dataToAdd.data_Y
            self.data_Z += dataToAdd.data_Z
        elif type(dataToAdd) == float:
            self.data_X += dataToAdd
            self.data_Y += dataToAdd
            self.data_Z += dataToAdd
        elif type(dataToAdd) == int:
            self.data_X += float(dataToAdd)
            self.data_Y += float(dataToAdd)
            self.data_Z += float(dataToAdd)

        return self

    def __sub__(self, dataToSub):
        """
        Des.
        """

        newData = InterfaceData(self.nPoint, self.mpiComm)
        if type(dataToSub) == type(self):
            newData.data_X = self.data_X - dataToSub.data_X
            newData.data_Y = self.data_Y - dataToSub.data_Y
            newData.data_Z = self.data_Z - dataToSub.data_Z
        elif type(dataToSub) == float:
            newData.data_X = self.data_X - dataToSub
            newData.data_Y = self.data_Y - dataToSub
            newData.data_Z = self.data_Z - dataToSub
        elif type(dataToSub) == int:
            newData.data_X = self.data_X - float(dataToSub)
            newData.data_Y = self.data_Y - float(dataToSub)
            newData.data_Z = self.data_Z - float(dataToSub)

        return newData

    def __rsub__(self, data):
        """
        des.
        """

        newData = -1*self + data

        return newData

    def __isub__(self, dataToSub):
        """
        Des.
        """
        if type(dataToSub) == type(self):
            self.data_X -= dataToSub.data_X
            self.data_Y -= dataToSub.data_Y
            self.data_Z -= dataToSub.data_Z
        elif type(dataToSub) == float:
            self.data_X -= dataToSub
            self.data_Y -= dataToSub
            self.data_Z -= dataToSub
        elif type(dataToSub) == int:
            self.data_X -= float(dataToSub)
            self.data_Y -= float(dataToSub)
            self.data_Z -= float(dataToSub)

        return self

    def __mul__(self, mulVal):
        """
        des.
        """

        newData = InterfaceData(self.nPoint, self.mpiComm)
        newData.data_X = self.data_X*mulVal
        newData.data_Y = self.data_Y*mulVal
        newData.data_Z = self.data_Z*mulVal

        return newData

    def __rmul__(self, mulVal):
        """
        des.
        """
        
        newData = self*mulVal

        return newData

    def __imul__(self, mulVal):
        """
        Des.
        """

        self.data_X *= mulVal
        self.data_Y *= mulVal
        self.data_Z *= mulVal

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
                 
#class InterfaceDisplacement(InterfaceData)
#class InterfaceLoad(InterfaceData)

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

    def __init__(self, FluidSolver, SolidSolver, nDim, computationType, mpiComm=None):
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

        self.solidInterfaceDisplacement = InterfaceData(self.ns, self.mpiComm)
        self.fluidInterfaceDisplacement = InterfaceData(self.nf, self.mpiComm)
        self.solidInterfaceLoads = InterfaceData(self.ns, self.mpiComm)
        self.fluidInterfaceLoads = InterfaceData(self.nf, self.mpiComm)

    #def createRTree(self, posX_array, posY_array, posZ_array):
        #"""
        #Description.
        #"""
        
        #from rtree import index
        #prop_index = index.Property()
        #prop_index.dimension = self.nDim
        #SpatialTree = index.Index(properties=prop_index)

        #for jVertex in range(self.ns):
        #    posX = posX_array[jVertex]
        #    posY = posY_array[jVertex]
        #    posZ = posZ_array[jVertex]
        #    if self.nDim == 2:
        #        SpatialTree.add(jVertex, (posX, posY))
        #    else:
        #        SpatialTree.add(jVertex, (posX, posY, posZ))

        #return SpatialTree

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
                self.solidInterfaceDisplacement[iGlobalVertex] = (localSolidInterfaceDisp_X[iVertex], localSolidInterfaceDisp_Y[iVertex], localSolidInterfaceDisp_Z[iVertex])

        self.solidInterfaceDisplacement.assemble()

    def getLoadsFromFluidSolver(self):
        """
        Des.
        """

        if self.myid in self.manager.getFluidInterfaceProcessors():
            localFluidInterfaceLoad_X, localFluidInterfaceLoad_Y, localFluidInterfaceLoad_Z = self.FluidSolver.getNodalLoads()
            for iVertex in range(self.nf_loc):
                iGlobalVertex = self.manager.getGlobalIndex('fluid', self.myid, iVertex)
                self.fluidInterfaceLoads[iGlobalVertex] = (localFluidInterfaceLoad_X[iVertex], localFluidInterfaceLoad_Y[iVertex], localFluidInterfaceLoad_Z[iVertex])

        self.fluidInterfaceLoads.assemble()

    def setLoadsToSolidSolver(self, time):
        """
        Des.
        """

        self.checkTotalLoad()

        # --- Redistribute the interpolated solid loads according to the partitions that own the solid interface --- #

        if self.mpiComm != None:
            localSize = self.solidInterfaceLoads.getXArray().shape[0]
            solidLoads_array_X_recon = mpiGatherv(self.solidInterfaceLoads.getXArray(), localSize, self.ns, self.mpiComm, 0)
            solidLoads_array_Y_recon = mpiGatherv(self.solidInterfaceLoads.getYArray(), localSize, self.ns, self.mpiComm, 0)
            solidLoads_array_Z_recon = mpiGatherv(self.solidInterfaceLoads.getZArray(), localSize, self.ns, self.mpiComm, 0)

            if self.myid == 0:
                for iProc in self.manager.getSolidInterfaceProcessors():
                    solidPhysicalInterfaceNodesDistribution = self.manager.getSolidPhysicalInterfaceNodesDistribution()
                    solidGlobalIndexRange = self.manager.getSolidGlobalIndexRange()
                    sendBuff_X = np.zeros(solidPhysicalInterfaceNodesDistribution[iProc])
                    sendBuff_Y = np.zeros(solidPhysicalInterfaceNodesDistribution[iProc])
                    sendBuff_Z = np.zeros(solidPhysicalInterfaceNodesDistribution[iProc])
                    globalIndex = solidGlobalIndexRange[iProc][iProc][0]
                    for iVertex in range(solidPhysicalInterfaceNodesDistribution[iProc]):
                        sendBuff_X[iVertex] = solidLoads_array_X_recon[globalIndex]
                        sendBuff_Y[iVertex] = solidLoads_array_Y_recon[globalIndex]
                        sendBuff_Z[iVertex] = solidLoads_array_Z_recon[globalIndex]
                        globalIndex += 1
                    self.mpiComm.Send(sendBuff_X, dest=iProc, tag = 1)
                    self.mpiComm.Send(sendBuff_Y, dest=iProc, tag = 2)
                    self.mpiComm.Send(sendBuff_Z, dest=iProc, tag = 3)
            if self.myid in self.manager.getSolidInterfaceProcessors():
                localSolidLoads_array_X = np.zeros(self.ns_loc)
                localSolidLoads_array_Y = np.zeros(self.ns_loc)
                localSolidLoads_array_Z = np.zeros(self.ns_loc)
                self.mpiComm.Recv(localSolidLoads_array_X, source=0, tag = 1)
                self.mpiComm.Recv(localSolidLoads_array_Y, source=0, tag = 2)
                self.mpiComm.Recv(localSolidLoads_array_Z, source=0, tag = 3)
                self.SolidSolver.applyNodalLoads(localSolidLoads_array_X, localSolidLoads_array_Y, localSolidLoads_array_Z, time)
        else:
            self.SolidSolver.applyNodalLoads(self.solidInterfaceLoads.getXArray(), self.solidInterfaceLoads.getYArray(), self.solidInterfaceLoads.getZArray(), time)

    def setDisplacementToFluidSolver(self, time):
        """
        Des.
        """

        self.checkConservation()

        if self.mpiComm != None:
            localSize = self.fluidInterfaceDisplacement.getXArray().shape[0]
            fluidInterface_array_DispX_recon = mpiGatherv(self.fluidInterfaceDisplacement.getXArray(), localSize, self.nf, self.mpiComm, 0)
            fluidInterface_array_DispY_recon = mpiGatherv(self.fluidInterfaceDisplacement.getYArray(), localSize, self.nf, self.mpiComm, 0)
            fluidInterface_array_DispZ_recon = mpiGatherv(self.fluidInterfaceDisplacement.getZArray(), localSize, self.nf, self.mpiComm, 0)
            haloNodesDisplacements = {}
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
                        sendBuff_X[iVertex] = fluidInterface_array_DispX_recon[globalIndex]
                        sendBuff_Y[iVertex] = fluidInterface_array_DispY_recon[globalIndex]
                        sendBuff_Z[iVertex] = fluidInterface_array_DispZ_recon[globalIndex]
                        globalIndex += 1
                    fluidHaloNodesList = self.manager.getFluidHaloNodesList()
                    fluidIndexing = self.manager.getFluidIndexing()
                    for key in fluidHaloNodesList[iProc].keys():
                        globalIndex = fluidIndexing[key]
                        sendBuffHalo[key] = (fluidInterface_array_DispX_recon[globalIndex], fluidInterface_array_DispY_recon[globalIndex], fluidInterface_array_DispZ_recon[globalIndex])
                    self.mpiComm.Send(sendBuff_X, dest=iProc, tag = 1)
                    self.mpiComm.Send(sendBuff_Y, dest=iProc, tag = 2)
                    self.mpiComm.Send(sendBuff_Z, dest=iProc, tag = 3)
                    self.mpiComm.send(sendBuffHalo, dest = iProc, tag=4)
            if self.myid in self.manager.getFluidInterfaceProcessors():
                localFluidInterface_array_DispX = np.zeros(self.nf_loc)
                localFluidInterface_array_DispY = np.zeros(self.nf_loc)
                localFluidInterface_array_DispZ = np.zeros(self.nf_loc)
                self.mpiComm.Recv(localFluidInterface_array_DispX, source=0, tag=1)
                self.mpiComm.Recv(localFluidInterface_array_DispY, source=0, tag=2)
                self.mpiComm.Recv(localFluidInterface_array_DispZ, source=0, tag=3)
                haloNodesDisplacements = self.mpiComm.recv(source=0, tag=4)
                self.FluidSolver.applyNodalDisplacements(localFluidInterface_array_DispX, localFluidInterface_array_DispY, localFluidInterface_array_DispZ, localFluidInterface_array_DispX, localFluidInterface_array_DispY, localFluidInterface_array_DispZ, haloNodesDisplacements, time)
        else:
            self.FluidSolver.applyNodalDisplacements(self.fluidInterfaceDisplacement.getXArray(), self.fluidInterfaceDisplacement.getYArray(), self.fluidInterfaceDisplacement.getZArray(), self.fluidInterfaceDisplacement.getXArray(), self.fluidInterfaceDisplacement.getYArray(), self.fluidInterfaceDisplacement.getZArray(), {}, time)

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

class MatchingMeshesInterpolator(InterfaceInterpolator):
    """
    Description.
    """

    def __init__(self, Manager, FluidSolver, SolidSolver, mpiComm = None):
        """
        Description
        """
        
        InterfaceInterpolator.__init__(self, Manager, FluidSolver, SolidSolver, mpiComm)

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

    def fillMatrix(self, solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc):
        """
        Des.
        """
            
        KDTree = self.createKDTree(solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z)
        #SolidSpatialTree = self.createRTree(solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z)
        localFluidInterface_array_X_init, localFluidInterface_array_Y_init, localFluidInterface_array_Z_init = self.FluidSolver.getNodalInitialPositions()
        for iVertex in range(self.nf_loc):
            posX = localFluidInterface_array_X_init[iVertex]
            posY = localFluidInterface_array_Y_init[iVertex]
            posZ = localFluidInterface_array_Z_init[iVertex]
            fluidPoint = np.array([posX, posY, posZ])
            #if self.nDim == 2:
            #    neighboors = list(SolidSpatialTree.nearest((posX, posY),1))
            #elif self.nDim == 3:
            #    neighboors = list(SolidSpatialTree.nearest((posX, posY, posZ),1))
            #jVertex = neighboors[0]
            distance, jVertex = KDTree.query(fluidPoint, 1)
            #NodeA = np.array([posX, posY, posZ])
            #NodeB = np.array([solidInterfaceBuffRcv_X[jVertex], solidInterfaceBuffRcv_Y[jVertex], solidInterfaceBuffRcv_Z[jVertex]])
            #distance = self.distance(NodeA, NodeB)
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

        self.H_T.mult(self.fluidInterfaceLoads.getDataX(), self.solidInterfaceLoads.getDataX())
        self.H_T.mult(self.fluidInterfaceLoads.getDataY(), self.solidInterfaceLoads.getDataY())
        self.H_T.mult(self.fluidInterfaceLoads.getDataZ(), self.solidInterfaceLoads.getDataZ())

    def interpolateSolidDisplacementOnFluidMesh(self):
        """
        Description.
        """

        self.H.mult(self.solidInterfaceDisplacement.getDataX(), self.fluidInterfaceDisplacement.getDataX())
        self.H.mult(self.solidInterfaceDisplacement.getDataY(), self.fluidInterfaceDisplacement.getDataY())
        self.H.mult(self.solidInterfaceDisplacement.getDataZ(), self.fluidInterfaceDisplacement.getDataZ())
    
#class NNInterpolator(InterfaceInterpolator):
    #"""
    #Description
    #"""

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
        self.writeInBGS = False
        
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
        
        self.solidInterfaceVelocity = InterfaceData(ns, self.mpiComm)
        self.solidInterfaceVelocitynM1 = InterfaceData(ns, self.mpiComm)
        self.solidInterfaceResidual = InterfaceData(ns, self.mpiComm)
    
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
        
        # --- Get the predicted (computed) solid interface displacement from the solid solver --- #
        predictedDisplacement = InterfaceData(ns, self.mpiComm)
        
        if self.myid in self.manager.getSolidInterfaceProcessors():
            localSolidInterfaceDisp_X, localSolidInterfaceDisp_Y, localSolidInterfaceDisp_Z = self.SolidSolver.getNodalDisplacements()
            for iVertex in range(self.manager.getNumberOfLocalSolidInterfaceNodes()):
                iGlobalVertex = self.manager.getGlobalIndex('solid', self.myid, iVertex)
                predictedDisplacement[iGlobalVertex] = (localSolidInterfaceDisp_X[iVertex], localSolidInterfaceDisp_Y[iVertex], localSolidInterfaceDisp_Z[iVertex])
        
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
                self.solidInterfaceVelocity[iGlobalVertex] = (localSolidInterfaceVel_X[iVertex], localSolidInterfaceVel_Y[iVertex], localSolidInterfaceVel_Z[iVertex])
                self.solidInterfaceVelocitynM1[iGlobalVertex] = (localSolidInterfaceVelNm1_X[iVertex], localSolidInterfaceVelNm1_Y[iVertex], localSolidInterfaceVelNm1_Z[iVertex])
        
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
    
    def bgsCoupling(self, timeIter, time):
        """
        Des
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
                self.FSIConv =    self.criterion.isVerified(self.errValue)
                
                # --- Relaxe the solid position --- #
                mpiPrint('\nProcessing interface displacements...\n', self.mpiComm)
                self.relaxSolidPosition()
            
            if self.writeInBGS == True:
                self.writeRealTimeData(timeIter, time, self.errValue)
            
            self.FSIIter += 1
        
        # --- Update the FSI history file --- #
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
            self.writeInBGS = True
            self.bgsCoupling(timeIter, time)
        
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
            self.bgsCoupling(timeIter, time)
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

class AlgortihmBGSAitkenRelax(AlgortihmBGSStaticRelax):
    
    def __init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, Criterion, nbFSIIterMax, deltaT, totTime, timeIterTreshold=-1, omegaMax=1.0, mpiComm=None):
        """
        Des.
        """
        
        AlgortihmBGSStaticRelax.__init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, Criterion, nbFSIIterMax, deltaT, totTime, timeIterTreshold, omegaMax, mpiComm)
        
        ns = self.interfaceInterpolator.getNs()
        self.solidInterfaceResidualnM1 = InterfaceData(ns, self.mpiComm)
    
    def setOmega(self):
        """
        Des.
        """
        
        if self.FSIIter != 0:
            # --- Compute the dynamic Aitken coefficient --- #
            deltaInterfaceResidual = self.solidInterfaceResidual - self.solidInterfaceResidualnM1
            
            prodScalRes_X, prodScalRes_Y, prodScalRes_Z = deltaInterfaceResidual.dot(self.solidInterfaceResidualnM1)
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
        self.solidInterfaceResidualnM1 = self.solidInterfaceResidual.copy()


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
