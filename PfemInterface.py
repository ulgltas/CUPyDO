#!/usr/bin/env python
# -*- coding: latin-1; -*-

# PfemInterface.py
# Python interface between the wrapper of PFEM solver and CUPyDO.
# Authors M.L. CERQUAGLIA
#
# COPYRIGHT (C) University of Liège, 2017.

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import os, os.path, sys, time, string

import math
import numpy as np
from FSICoupler import FluidSolver

# ----------------------------------------------------------------------
#  PfemSolver class
# ----------------------------------------------------------------------

class PfemSolver(FluidSolver):
    def __init__(self, testname, bndno, dt):
        
        print '\n***************************** Initializing Pfem *****************************'
        
        self.testname = testname  # string (name of the module of the fluid model)
        self.bndno = bndno        # phygroup# of the f/s interface
        
        # internal vars
        self.vnods = []           # dict of interface nodes / prescribed velocities
        self.t1      = 0.0        # last reference time        
        self.t2      = 0.0        # last calculated time
                 
        # loads the python module
        #load(self.testname)         # use toolbox.utilities
        exec("import %s" % self.testname)
        exec("module = %s" % self.testname) # link to Pfem object

        # create an instance of Pfem
        self.pfem = module.getPfem()
        self.realTimeExtractorsList = module.getRealTimeExtractorsList(self.pfem)
        
        # retrieve the f/s boundary and the related nodes
        gr = self.pfem.w.Group(self.pfem.msh, bndno)
        
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
        self.V = self.pfem.w.DoubleVector()
        self.V0 = self.pfem.w.DoubleVector()
        self.u = self.pfem.w.DoubleVector()
        self.v = self.pfem.w.DoubleVector()
        self.p = self.pfem.w.DoubleVector()
        self.velocity = self.pfem.w.DoubleVector()
        
        self.pfem.scheme.t = 0.
        self.pfem.scheme.dt = dt
        self.pfem.scheme.init(self.V,self.V0,self.u,self.v,self.p,self.velocity)
        
        self.runOK = True
        
        FluidSolver.__init__(self)
        
        self.displ_x_Nm1 = np.zeros((self.nPhysicalNodes))
        self.displ_y_Nm1 = np.zeros((self.nPhysicalNodes))
        self.displ_z_Nm1 = np.zeros((self.nPhysicalNodes))
        
    def run(self, t1, t2):
        """
        calculates one increment from t1 to t2.
        """
        
        self.pfem.scheme.setNextStep(self.V,self.V0,self.u,self.v,self.p,self.velocity)
        self.runOK = self.pfem.scheme.runOneTimeStep(self.V,self.V0,self.u,self.v,self.p,self.velocity)
        
        self.__setCurrentState()
    
    def getNodalInitialPositions(self):
        
        x0 = np.zeros(len(self.vnods))
        y0 = np.zeros(len(self.vnods))
        z0 = np.zeros(len(self.vnods))
        
        for i in range(len(self.vnods)):
            node = self.vnods[i]                 
            x0[i] = node.pos0.x[0]
            y0[i] = node.pos0.x[1]
            z0[i] = 0.
        
        return x0, y0, z0
    
    def __setCurrentState(self):
        
        fx = np.zeros(len(self.vnods))
        fy = np.zeros(len(self.vnods))
        fz = np.zeros(len(self.vnods))
        
        for i in range(len(self.vnods)):
            node = self.vnods[i]
            fx[i] = -node.Fint.x[0]
            fy[i] = -node.Fint.x[1]
            fz[i] = 0.
        
        self.nodalLoad_X = fx
        self.nodalLoad_Y = fy
        self.nodalLoad_Z = fz
    
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
        for no in self.vnods.iterkeys():
            node = self.vnods[no]                 
            vx = -0.5 # current vx          
            vy = 0. # current vy
            px = node.pos0.x[0] + vx*self.pfem.scheme.dt # current x         
            py = node.pos0.x[1] + vy*self.pfem.scheme.dt# current y              
            out[no] = (px,py,vx,vy)
            
        self.applydefo(out)
    
    def applyNodalDisplacements(self, dx, dy, dz, dx_nM1, dy_nM1, dz_nM1, haloNodesDisplacements,time):
        """
        prescribes given nodal positions and velocities coming from solid solver to node #no
        """
        if self.pfem.scheme.t < time:
            self.pfem.scheme.resetNodalPositions()
        
        for i in range(len(self.vnods)):
            node = self.vnods[i]                 
            node.imposedU = (dx[i] - self.displ_x_Nm1[i])/self.pfem.scheme.dt
            node.imposedV = (dy[i] - self.displ_y_Nm1[i])/self.pfem.scheme.dt
        
    def update(self, dt):
        self.pfem.scheme.t+=dt
        self.pfem.scheme.nt+=1
        self.pfem.scheme.updateSolutionVectors(self.V,self.V0,self.u,self.v,self.p,self.velocity)
        
        
        #---
        for i in range(len(self.vnods)):
            displ_x = self.displ_x_Nm1[i] + self.u[self.vnods[i].rowU]*dt
            displ_y = self.displ_y_Nm1[i] + self.v[self.vnods[i].rowU]*dt
            self.displ_x_Nm1[i] = displ_x
            self.displ_y_Nm1[i] = displ_y
        # ---
        
    def save(self, nt):
        if nt%self.pfem.scheme.savefreq==0:
            self.pfem.scheme.archive()
        if not self.pfem.gui==None:
            self.pfem.scheme.vizu(self.u,self.v,self.p)
        
    def initRealTimeData(self):
        """
        Des.
        """
        
        for extractor in self.realTimeExtractorsList:
            extractorName = extractor.dofName
            solFile = open(extractorName + '.ascii', "w")
            solFile.write("Time\tnIter\tValue\n")
            solFile.close() #Should we keep it open?
    
    def saveRealTimeData(self, time, nFSIIter):
        """
        Des.
        """
        
        for extractor in self.realTimeExtractorsList:
            data = extractor.extract()
            extractorName = extractor.dofName
            solFile = open(extractorName + '.ascii', "a")
            buff = str()
            for ii in range(data.size()):
                buff = buff + '\t' + str(data[ii])
            solFile.write(str(time) + '\t' + str(nFSIIter) + buff + '\n')
            solFile.close()
    
    def printRealTimeData(self, time, nFSIIter):
        """
        Des.
        """
        
        for extractor in self.realTimeExtractorsList:
            data = extractor.extract()
            extractorName = extractor.dofName
            buff = str()
            for ii in range(data.size()):
                buff = buff + '\t' + str(data[ii])
            toPrint = 'RES-FSI-' + extractorName + ': ' + buff + '\n'
            print toPrint
    
    def remeshing(self):
        self.pfem.scheme.remeshing(self.V,self.V0,self.p)
        self.pfem.scheme.updateData()
    
    def exit(self):
        """
        Exits the Pfem solver.
        """
        if not self.pfem.gui==None:
            self.pfem.gui.save2vtk()
        self.pfem.scheme.exit()
        
        print("***************************** Exit Pfem *****************************")
