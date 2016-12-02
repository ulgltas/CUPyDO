#!/usr/bin/env python
# -*- coding: latin-1; -*-

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import math
from toolbox.utilities import *
import toolbox.fac as fac
from wrap import *
import numpy as np
from fsi import SolidSolver

# ----------------------------------------------------------------------
#  Nodal Load class
# ----------------------------------------------------------------------

class NLoad:
    """
    Nodal load
    """
    def __init__(self, val1, t1, val2, t2):
        self.val1 = val1
        self.t1 = t1
        self.val2 = val2
        self.t2 = t2
    def __call__(self, time):
        return self.val1 + (time-self.t1)/(self.t2-self.t1)*(self.val2-self.val1)
    def nextstep(self):
        self.t1 = self.t2
        self.val1 = self.val2

# ----------------------------------------------------------------------
#  MtfSolver class
# ----------------------------------------------------------------------
               
class MtfSolver(SolidSolver):
    def __init__(self, testname):
        """
        des.
        """

        # --- Load the Python module --- #
        self.testname = testname            # string (name of the module of the solid model)
	#load(self.testname)                # loads the python module and creates mtf/workspace
	exec("import %s" % self.testname)
        exec("module = %s" % self.testname)

        # --- Create an instance of Metafor --- #
        self.metafor = None                   # link to Metafor objet
        p = {}                                # parameters (todo)
        #self.metafor = instance(p)           # creates an instance of the metafor model
	self.metafor = module.getMetafor(p)
        self.realTimeExtractorsList = module.getRealTimeExtractorsList(self.metafor)
	p = module.params(p)

        # --- Internal variables --- #
	self.neverRun = True            # bool True until the first Metafor run is completed then False
        self.fnods = {}                 # dict of interface nodes / prescribed forces
        self.t1      = 0.0              # last reference time        
        self.t2      = 0.0              # last calculated time
        self.nbFacs = 0                 # number of existing Facs
        self.saveAllFacs = True         # True: the Fac corresponding to the end of the time step is conserved, False: Facs are erased at the end of each time step
        self.runOK = True

        # --- Retrieves the f/s boundary and the related nodes --- #
	self.bndno = p['bndno']                                                 # physical group of the f/s interface
        #self.extractNode = p['extractNode']                                     # physical group of the node to be extracted (output)
        self.groupset = self.metafor.getDomain().getGeometry().getGroupSet()
        self.gr = self.groupset(self.bndno)
        self.nLocalInterfaceNodes = self.gr.getNumberOfMeshPoints()
        self.nLocalInterfaceHaloNode = 0
        self.nLocalInterfacePhysicalNodes = self.gr.getNumberOfMeshPoints()                     # number of node at the f/s boundary

        SolidSolver.__init__(self)

        # --- Builds a list (dict) of interface nodes and creates the nodal prescribed loads --- #
        loadingset = self.metafor.getDomain().getLoadingSet()
        for i in range(self.nLocalInterfacePhysicalNodes):
            node = self.gr.getMeshPoint(i)
            no = node.getNo()  
            fx = NLoad(self.t1, 0.0, self.t2, 0.0)
            fy = NLoad(self.t1, 0.0, self.t2, 0.0)
	    fz = NLoad(self.t1, 0.0, self.t2, 0.0)
            self.fnods[no] = (node, fx, fy, fz)
            fctx = PythonOneParameterFunction(fx)
            fcty = PythonOneParameterFunction(fy)
	    fctz = PythonOneParameterFunction(fz)
            loadingset.define(node, Field1D(TX,GF1), 1.0, fctx) 
            loadingset.define(node, Field1D(TY,GF1), 1.0, fcty)
	    loadingset.define(node, Field1D(TZ,GF1), 1.0, fctz)

        # --- Create the array for external communication (displacement, velocity and velocity at the previous time step --- #
        self.__setCurrentState()
        self.interfaceVel_XNm1 = self.interfaceVel_X.copy()
        self.interfaceVel_YNm1 = self.interfaceVel_Y.copy()
        self.interfaceVel_ZNm1 = self.interfaceVel_Z.copy()

    def run(self, t1, t2):
        """
        calculates one increment from t1 to t2.
        """
        if(self.neverRun):
            self.__firstRun(t1, t2)
            self.neverRun=False
        else:
            self.__nextRun(t1, t2)
        self.t1 = t1
        self.t2 = t2

        self.__setCurrentState()

    def __firstRun(self, t1, t2):
        """
        performs a first run of metafor with all the required preprocessing.
        """
        # this is the first run - initialize the timestep manager of metafor
        tsm = self.metafor.getTimeStepManager()
        dt    = t2-t1  # time-step size
        dt0   = dt     # initial time step
        dtmax = dt     # maximum size of the time step
        tsm.setInitialTime(t1, dt0)
        tsm.setNextTime(t2, 1, dtmax)
        # launches metafor from t1 to t2
        #meta()                  # use toolbox.utilities
        log = LogFile("resFile.txt")
        self.runOK = self.metafor.getTimeIntegration().integration()
        # at this stage, 2 archive files have been created in the workspace

    def __nextRun(self, t1, t2):
        """
        performs one time increment (from t1 to t2) of the solid model.
        this increment is a full metafor run and it may require more than 1 time step.
        """
        if self.t1==t1:
            # rerun from t1
            if self.t2!=t2:
                raise Exception("bad t2 (%f!=%f)" % (t2, self.t2)) 
            
            loader = fac.FacManager(self.metafor)
            nt = loader.lookForFile(self.nbFacs) #(0)
            loader.eraseAllFrom(nt)
            self.runOK = self.metafor.getTimeIntegration().restart(nt)
        else:
            # new time step
            tsm = self.metafor.getTimeStepManager()
            dt=t2-t1
            dtmax=dt
            tsm.setNextTime(t2, 1, dtmax/2)    # forces at least 2 time increments   
            
            loader = fac.FacManager(self.metafor)
            print 'self.nbFacs = ', self.nbFacs
            nt1 = loader.lookForFile(self.nbFacs) #(0)
            nt2 = loader.lookForFile(self.nbFacs+1) #(1)
            if not self.saveAllFacs:
                loader.erase(nt1) # delete first fac
            self.runOK = self.metafor.getTimeIntegration().restart(nt2)
            if self.saveAllFacs:
                self.nbFacs+=1

    def __setCurrentState(self):
        """
        Des.
        """

        for ii in range(self.nLocalInterfacePhysicalNodes):
            node = self.gr.getMeshPoint(ii)
            posX = node.getPos(Configuration().currentConf).get1() # current x coord
            posY = node.getPos(Configuration().currentConf).get2() # current y coord
            posZ = node.getPos(Configuration().currentConf).get3() # current z coord
            velX = node.getValue(Field1D(TX,GV))			#current vel x
	    velY = node.getValue(Field1D(TY,GV))			#current vel y
	    velZ = node.getValue(Field1D(TZ,GV))                  #current vel z
            self.interfaceDisp_X[ii] = node.getValue(Field1D(TX,RE))
            self.interfaceDisp_Y[ii] = node.getValue(Field1D(TY,RE))
	    self.interfaceDisp_Z[ii] = node.getValue(Field1D(TZ,RE))
            self.interfaceVel_X[ii] = velX
            self.interfaceVel_Y[ii] = velY
            self.interfaceVel_Z[ii] = velZ

    def getInterfaceNodalInitialPositions(self):
        """
        Description.
        """

        interfaceInitialPos_X = np.zeros(self.nLocalInterfacePhysicalNodes)
        interfaceInitialPos_Y = np.zeros(self.nLocalInterfacePhysicalNodes)
        interfaceInitialPos_Z = np.zeros(self.nLocalInterfacePhysicalNodes)

        for ii in range(self.nLocalInterfacePhysicalNodes):
            node = self.gr.getMeshPoint(ii)
            interfaceInitialPos_X[ii] = node.getPos0().get1()   
            interfaceInitialPos_Y[ii] = node.getPos0().get2()
            interfaceInitialPos_Z[ii] = node.getPos0().get3() 

        return (interfaceInitialPos_X, interfaceInitialPos_Y, interfaceInitialPos_Z)

    def getInterfaceNodeIndex(self, iVertex):
	"""
	Returns the index (identifier) of the iVertex^th interface node.
	"""
	node = self.gr.getMeshPoint(iVertex)
        no = node.getNo()

	return no

    def fakeFluidSolver(self, time):
        """
        calculate some dummy loads as a function of timestep.
        these loads should be replaced by the fluid solver in practice.
        for each node, the fsi solver may call the "solid.applyload" function.
        """
        # calculate L (max length along x)
        xmin=1e10
        xmax=-1e10
        for no in self.fnods.iterkeys():
            node,fx,fy = self.fnods[no]
            px = node.getPos0().get1()
            if px<xmin: xmin=px
            if px>xmax: xmax=px
        #print "(xmin,xmax)=(%f,%f)" % (xmin,xmax)
        L = xmax-xmin
        #print "L=%f" % L
    
        # loop over node#
        for no in self.fnods.iterkeys():
            node,fx,fy = self.fnods[no]  # retreive data of node #no
            px = node.getPos0().get1()     # x coordinate of the node       
            valx = 0.0 
            valy = -3e-4*time*math.sin(8*math.pi*px/L) # dummy fct
            self.applyload(no, valx, valy, time)

    def applyLoad(self, load_X, load_Y, load_Z, time):
        """
        Des.
        """

        for ii in range(self.nLocalInterfacePhysicalNodes):
            node = self.gr.getMeshPoint(ii)
            no = node.getNo()
            node,fx,fy,fz = self.fnods[no]
            fx.val2 = load_X[ii]
            fy.val2 = load_Y[ii]
            fz.val2 = load_Z[ii]
            fx.t2 = time
            fy.t2 = time
            fz.t2 = time

    def __nextStep(self):
        """
        Des.
        """

        for no in self.fnods.iterkeys():
            node, fx, fy, fz = self.fnods[no]
            fx.nextstep()
            fy.nextstep()
            fz.nextstep()

    def update(self):
	"""
	Pushes back the current state in the past (previous state) before going to the next time step.
	"""

        SolidSolver.update(self)
        
        self.__nextStep()

    def initRealTimeData(self):
        """
        Des.
        """

        for extractor in self.realTimeExtractorsList:
            #extractorName = extractor.name()
            extractorName = 'TestExtractor'
            solFile = open(extractorName + '.ascii', "w")
            solFile.write("Time\tIteration no.\tValue\n")
            solFile.close()

    def saveRealTimeData(self, time, nFSIIter):
        """
        Des.
        """

        for extractor in self.realTimeExtractorsList:
            data = extractor.extract()
            #extractorName = extractor.name()
            extractorName = 'TestExtractor'
            solFile = open(extractorName + '.ascii', "a")
            buff = str()
            for ii in range(data.size()):
              buff = buff + '\t' + str(data[ii])
            solFile.write(str(time) + '\t' + str(nFSIIter) + buff + '\n')
            solFile.close()

    def exit(self):
	"""
	Exits the Metafor solver.
	"""
	print("***************************** Exit Metafor *****************************")
