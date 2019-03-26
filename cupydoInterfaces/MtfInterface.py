#! /usr/bin/env python
# -*- coding: latin-1; -*-

''' 

Copyright 2018 University of Li�ge

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 

MtfInterface.py
Python interface between the wrapper of Metafor and CUPyDO.
Authors R. BOMAN, M.L. CERQUAGLIA, D. THOMAS

'''

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import os, os.path, sys, time, string

import math
from toolbox.utilities import *
import toolbox.fac as fac
from wrap import *
import numpy as np
from cupydo.genericSolvers import SolidSolver

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
        theValue = self.val1 + (time-self.t1)/(self.t2-self.t1)*(self.val2-self.val1)
        return theValue
    def nextstep(self):
        self.t1 = self.t2
        self.val1 = self.val2

# ----------------------------------------------------------------------
#  MtfSolver class
# ----------------------------------------------------------------------
               
class MtfSolver(SolidSolver):
    def __init__(self, testname, computationType):
        """
        des.
        """
        
        print '\n***************************** Initializing Metafor *****************************'
        
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
        self.Tnods = {}
        self.t1      = 0.0              # last reference time        
        self.t2      = 0.0              # last calculated time
        self.nbFacs = 0                 # number of existing Facs
        self.saveAllFacs = False         # True: the Fac corresponding to the end of the time step is conserved, False: Facs are erased at the end of each time step
        self.runOK = True
        self.computationType = computationType  # computation type : steady (default) or unsteady

        # --- Retrieves the f/s boundary and the related nodes --- #
        self.bndno = p['bndno']                                                 # physical group of the f/s interface
        self.groupset = self.metafor.getDomain().getGeometry().getGroupSet()
        self.gr = self.groupset(self.bndno)
        self.nNodes = self.gr.getNumberOfMeshPoints()
        self.nHaloNode = 0
        self.nPhysicalNodes = self.gr.getNumberOfMeshPoints()                     # number of node at the f/s boundary


        # --- Builds a list (dict) of interface nodes and creates the nodal prescribed loads --- #
        loadingset = self.metafor.getDomain().getLoadingSet()
        for i in range(self.nPhysicalNodes):
            node = self.gr.getMeshPoint(i)
            no = node.getNo()
            fx = NLoad(self.t1, 0.0, self.t2, 0.0)
            fy = NLoad(self.t1, 0.0, self.t2, 0.0)
            fz = NLoad(self.t1, 0.0, self.t2, 0.0)
            Temp = NLoad(self.t1, 0.0, self.t2, 0.0)
            self.fnods[no] = (node, fx, fy, fz)
            self.Tnods[no] = (node, Temp)
            fctx = PythonOneParameterFunction(fx)
            fcty = PythonOneParameterFunction(fy)
            fctz = PythonOneParameterFunction(fz)
            fctTemp = PythonOneParameterFunction(Temp)
            loadingset.define(node, Field1D(TX,GF1), 1.0, fctx) 
            loadingset.define(node, Field1D(TY,GF1), 1.0, fcty)
            loadingset.define(node, Field1D(TZ,GF1), 1.0, fctz)
            #loadingset.define(node, Field1D(TO, RE), 1.0, fctTemp, 0)
            #loadingset.define(node, Field1D(TO, AB), 1.0, fctTemp)
        
        # --- Create the array for external communication (displacement, velocity and velocity at the previous time step) --- #
        
        SolidSolver.__init__(self)
        
        self.__setCurrentState(True)
        self.nodalVel_XNm1 = self.nodalVel_X.copy()
        self.nodalVel_YNm1 = self.nodalVel_Y.copy()
        self.nodalVel_ZNm1 = self.nodalVel_Z.copy()
        
        # Last build operation
        self.metafor.getDomain().build() # NB: necessary to complete Metafor initialization!
        
        self.initRealTimeData() #NB: to be called after self.metafor.getDomain().build() otherwise no proper extractors usage!
        
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

        self.__setCurrentState(False)

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
            tsm.setNextTime(t2, 1, dtmax)  
            
            loader = fac.FacManager(self.metafor)
            nt1 = loader.lookForFile(self.nbFacs) #(0)
            nt2 = loader.lookForFile(self.nbFacs+1) #(1)
            if not self.saveAllFacs:
                loader.erase(nt1) # delete first fac
            self.runOK = self.metafor.getTimeIntegration().restart(nt2)
            if self.saveAllFacs:
                self.nbFacs+=1

    def __setCurrentState(self, val_init):
        """
        Des.
        """

        for ii in range(self.nPhysicalNodes):
            node = self.gr.getMeshPoint(ii)
            posX = node.getPos(Configuration().currentConf).get1()    # current x coord
            posY = node.getPos(Configuration().currentConf).get2()    # current y coord
            posZ = node.getPos(Configuration().currentConf).get3()    # current z coord
            velX = node.getValue(Field1D(TX,GV))                      # current vel x
            velY = node.getValue(Field1D(TY,GV))                      # current vel y
            velZ = node.getValue(Field1D(TZ,GV))                      # current vel z
            self.nodalDisp_X[ii] = node.getValue(Field1D(TX,RE))  # current x displacement
            self.nodalDisp_Y[ii] = node.getValue(Field1D(TY,RE))  # current y displacement
            self.nodalDisp_Z[ii] = node.getValue(Field1D(TZ,RE))  # current z displacement
            self.nodalVel_X[ii] = velX
            self.nodalVel_Y[ii] = velY
            self.nodalVel_Z[ii] = velZ
            #self.nodalTemperature[ii] = node.getValue(Field1D(TO,RE))
            #self.nodalTemperature[ii] = node.getValue(Field1D(TO,AB))
            if not val_init:
                HF_X_extractor = IFNodalValueExtractor(node, IF_FLUX_X)
                HF_Y_extractor = IFNodalValueExtractor(node, IF_FLUX_Y)
                HF_Z_extractor = IFNodalValueExtractor(node, IF_FLUX_Z)
                HF_X_data = HF_X_extractor.extract()
                HF_Y_data = HF_Y_extractor.extract()
                HF_Z_data = HF_Z_extractor.extract()
                self.nodalHeatFlux_X[ii] = HF_X_data[0]
                self.nodalHeatFlux_Y[ii] = HF_Y_data[0]
                self.nodalHeatFlux_Z[ii] = HF_Z_data[0]
                pos0X = node.getPos0().get1()
                pos0Y = node.getPos0().get2()
                pos0Z = node.getPos0().get3()

    def getNodalInitialPositions(self):
        """
        Description.
        """

        nodalInitialPos_X = np.zeros(self.nPhysicalNodes)
        nodalInitialPos_Y = np.zeros(self.nPhysicalNodes)
        nodalInitialPos_Z = np.zeros(self.nPhysicalNodes)

        for ii in range(self.nPhysicalNodes):
            node = self.gr.getMeshPoint(ii)
            nodalInitialPos_X[ii] = node.getPos0().get1()   
            nodalInitialPos_Y[ii] = node.getPos0().get2()
            nodalInitialPos_Z[ii] = node.getPos0().get3() 

        return (nodalInitialPos_X, nodalInitialPos_Y, nodalInitialPos_Z)

    def getNodalIndex(self, iVertex):
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
        
        valx = []
        valy = []
        valz = []
        
        # calculate L (max length along x)
        xmin=1e10
        xmax=-1e10
        
        for ii in range(self.nPhysicalNodes):
            node = self.gr.getMeshPoint(ii)
            no = node.getNo()
            node,fx,fy,fz = self.fnods[no]
            px = node.getPos0().get1()
            if px<xmin: xmin=px
            if px>xmax: xmax=px
        L = xmax-xmin
    
        # loop over node#
        for ii in range(self.nPhysicalNodes):
            node = self.gr.getMeshPoint(ii)
            no = node.getNo()
            node,fx,fy,fz = self.fnods[no]
            px = node.getPos0().get1()
            valx.append(0.)
            valy.append(-3e-4*time*math.sin(8*math.pi*px/L)) # dummy fct
            valz.append(0.)
        
        self.applyNodalLoads(valx, valy, valz, time)

    def applyNodalLoads(self, load_X, load_Y, load_Z, val_time):
        """
        Des.
        """

        for ii in range(self.nPhysicalNodes):
            node = self.gr.getMeshPoint(ii)
            no = node.getNo()
            node,fx,fy,fz = self.fnods[no]
            fx.val2 = load_X[ii]
            fy.val2 = load_Y[ii]
            fz.val2 = load_Z[ii]
            fx.t2 = val_time
            fy.t2 = val_time
            fz.t2 = val_time

    def applyNodalTemperatures(self, Temperature, val_time):
        """
        Des.
        """

        for ii in range(self.nPhysicalNodes):
            node = self.gr.getMeshPoint(ii)
            no = node.getNo()
            node, Temp = self.Tnods[no]
            Temp.val2 = Temperature[ii]
            Temp.t2 = val_time

    def __nextStep(self):
        """
        Des.
        """
        
        for no in self.fnods.iterkeys():
            node, fx, fy, fz = self.fnods[no]
            fx.nextstep()
            fy.nextstep()
            fz.nextstep()

        for no in self.Tnods.iterkeys():
            node, Temp = self.Tnods[no]
            Temp.nextstep()

    def update(self):
        """
        Pushes back the current state in the past (previous state) before going to the next time step.
        """

        SolidSolver.update(self)

        self.__nextStep()

    def bgsUpdate(self):
        """
        Des.
        """

        if self.computationType == 'steady':
            self.__nextStep()

    def initRealTimeData(self):
        """
        Des.
        """
        
        for extractor in self.realTimeExtractorsList:
            data = extractor.extract()
            extractorName = extractor.buildName()
            solFile = open(extractorName + '.ascii', "w")
            solFile.write("{0:>12s}   {1:>12s}".format("Time", "FSI_Iter"))
            for ii in range(len(data)):
                solFile.write("   {0:>12s}".format("Value_"+str(ii)))
            solFile.write("\n")
            solFile.close()

    def saveRealTimeData(self, time, nFSIIter):
        """
        Des.
        """
        
        for extractor in self.realTimeExtractorsList:
            data = extractor.extract()
            extractorName = extractor.buildName()
            solFile = open(extractorName + '.ascii', "a")
            solFile.write("{0:12.6f}   {1:12d}".format(time, nFSIIter))
            for d in data:
                solFile.write("   {0:12.6f}".format(d))
            solFile.write("\n")
            solFile.close()

    def printRealTimeData(self, time, nFSIIter):
        """
        Des.
        """
        
        for extractor in self.realTimeExtractorsList:
            data = extractor.extract()
            extractorName = extractor.buildName()
            buff = str()
            for d in data:
                buff = buff + '\t' + str(d)
            toPrint = 'RES-FSI-' + extractorName + ': ' + buff
            print toPrint
    
    def exit(self):
        """
        Exits the Metafor solver.
        """
        print("\n***************************** Exit Metafor *****************************")
