#! /usr/bin/env python3
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

Pfem.py
Python interface between the wrapper of PFEM solver and CUPyDO.
Authors M.L. CERQUAGLIA

'''

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import os, os.path, sys, time, string

import math
import numpy as np
from ..utilities import titlePrint
from ..genericSolvers import FluidSolver

# ----------------------------------------------------------------------
#  Pfem solver interface class
# ----------------------------------------------------------------------

class Pfem(FluidSolver):
    def __init__(self, testname, nthreads, dt):
        
        titlePrint('Initializing Pfem')
        
        self.testname = testname  # string (name of the module of the fluid model)
        
        # internal vars
        self.vnods = []           # dict of interface nodes / prescribed velocities
        self.t1      = 0.0        # last reference time        
        self.t2      = 0.0        # last calculated time
                 
        # loads the python module
        #load(self.testname)         # use toolbox.utilities
        exec("import %s" % self.testname, globals())
        exec("module = %s" % self.testname, globals()) # link to Pfem object

        # create an instance of Pfem
        self.pfem = module.getPfem()
        self.realTimeExtractorsList = module.getRealTimeExtractorsList(self.pfem)

        # Current mesh state VTK extractor
        import pfem.tools.link2vtk as v
        self.extractor = v.Link2VTK(self.pfem.msh, self.pfem.scheme, False, True)
        
        # retrieve the f/s boundary and the related nodes
        gr = self.pfem.w.Group(self.pfem.msh, self.pfem.bndno)
        
        # builds a list (dict) of interface nodes
        nods = {}
        for e in gr.tag.elems:
            for n in e.nodes:
                no = n.no
                nods[no] = n
        
        self.vnods = list(nods.values())
        
        self.nNodes = len(self.vnods)
        self.nHaloNode = 0    # numbers of nodes at the f/s interface (halo)
        self.nPhysicalNodes = self.nNodes - self.nHaloNode                        # numbers of nodes at the f/s interface (physical)
        
        # Pfem scheme initialization
        #self.u = self.pfem.w.DoubleVector()
        #self.v = self.pfem.w.DoubleVector()
        #self.p = self.pfem.w.DoubleVector()
        #self.velocity = self.pfem.w.DoubleVector()
        #self.matID = self.pfem.w.IntVector()
        
        self.pfem.scheme.t = 0.
        self.pfem.scheme.dt = dt
        self.pfem.scheme.init()
        # [AC] I moved the following 3 lines from the test cases definition to the interface
        self.pfem.scheme.nthreads = nthreads
        
        self.runOK = True
        
        FluidSolver.__init__(self)
        
        self.displ_x_Nm1 = np.zeros((self.nPhysicalNodes))
        self.displ_y_Nm1 = np.zeros((self.nPhysicalNodes))
        self.displ_z_Nm1 = np.zeros((self.nPhysicalNodes))
        
    def run(self, t1, t2):
        """
        calculates one increment from t1 to t2.
        """
        
        self.pfem.scheme.setNextStep()
        self.runOK = self.pfem.scheme.runOneTimeStep()
        
        self.__setCurrentState()
    
    def getNodalInitialPositions(self):
        
        x0 = np.zeros(len(self.vnods))
        y0 = np.zeros(len(self.vnods))
        z0 = np.zeros(len(self.vnods))
        
        for i in range(len(self.vnods)):
            node = self.vnods[i]                 
            x0[i] = node.posN.x[0]
            y0[i] = node.posN.x[1]
            z0[i] = 0.
        
        return x0, y0, z0

    def __setCurrentState(self):
        
        fx = np.zeros(len(self.vnods))
        fy = np.zeros(len(self.vnods))
        fz = np.zeros(len(self.vnods))

        for i in range(len(self.vnods)):
            node = self.vnods[i]
            fx[i] = -(node.fIne.x[0] + node.fInt.x[0] - node.fExt.x[0])
            fy[i] = -(node.fIne.x[1] + node.fInt.x[1] - node.fExt.x[1])
            fz[i] = 0.
        
        self.nodalLoad_X = fx
        self.nodalLoad_Y = fy
        self.nodalLoad_Z = fz

        if self.pfem.pbl.isAxisymmetric: # Rescale to 1 radian for Metafor (B-J)

            self.nodalLoad_X /= (2.0*np.pi)
            self.nodalLoad_Y /= (2.0*np.pi)
    
    def getNodalIndex(self, iVertex):
        """
        Returns the index (identifier) of the iVertex^th interface node.
        """
        node = self.vnods[iVertex]
        no = node.no
        
        return no
    
    def fakeSolidSolver(self, time):
        """
        calculate some dummy positions and velocities as a function of timestep.
        it should be replaced by the solid solver in practice.
        for each node, the fsi solver may call the "fluid.applyposition" function.
        """
        
        out = {}
        for no in self.vnods.keys():
            node = self.vnods[no]                 
            vx = -0.5 # current vx          
            vy = 0. # current vy
            px = node.posN.x[0] + vx*self.pfem.scheme.dt # current x         
            py = node.posN.x[1] + vy*self.pfem.scheme.dt# current y              
            out[no] = (px,py,vx,vy)
            
        self.applydefo(out)
    
    def applyNodalDisplacements(self, dx, dy, dz, dx_nM1, dy_nM1, dz_nM1, haloNodesDisplacements,time):
        """
        prescribes given nodal positions and velocities coming from solid solver to node #no
        """
        # if self.pfem.scheme.t < time:
        #     self.pfem.scheme.resetNodalPositions()
        
        for i in range(len(self.vnods)):
            node = self.vnods[i]                 
            node.imposedU = (dx[i] - self.displ_x_Nm1[i])/self.pfem.scheme.dt
            node.imposedV = (dy[i] - self.displ_y_Nm1[i])/self.pfem.scheme.dt
        
    def update(self, dt):
        self.pfem.scheme.t+=dt
        self.pfem.scheme.nt+=1
        self.pfem.scheme.updateSolutionVectors()
        #self.pfem.scheme.updateMatIDVector()
        #self.u = self.pfem.scheme.u
        #self.v = self.pfem.scheme.v
        
        #---
        for i in range(len(self.vnods)):
            displ_x = self.displ_x_Nm1[i] + self.pfem.scheme.u[self.vnods[i].rowU]*dt
            displ_y = self.displ_y_Nm1[i] + self.pfem.scheme.v[self.vnods[i].rowU]*dt
            self.displ_x_Nm1[i] = displ_x
            self.displ_y_Nm1[i] = displ_y
        # ---
        
    def save(self, nt):
            self.pfem.scheme.archive()
            self.pfem.scheme.vizu()
            self.extractor.saveVTK(nt,self.pfem.scheme.dt)
        
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
            extractorName = extractor.dofName
            buff = str()
            for d in data:
                buff = buff + '\t' + str(d)
            toPrint = 'RES-FSI-' + extractorName + ': ' + buff + '\n'
            print(toPrint)
    
    def remeshing(self):
        self.pfem.scheme.remeshing()
        self.pfem.scheme.updateNodalRefValues()
    
    def exit(self):
        """
        Exits the Pfem solver.
        """

        self.pfem.scheme.exit()
        titlePrint("Exit Pfem")
