#!/usr/bin/env python
# -*- coding: latin-1; -*-

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import os, os.path, sys, time, string

pfemPath = 'D:/PFEM/pfemB/bin/Release'
pfemToolsPath = 'D:/PFEM/pfem/tools' 
sys.path.append(pfemPath)
sys.path.append(pfemToolsPath)

import math
import numpy as np

# ----------------------------------------------------------------------
#  PfemSolver class
# ----------------------------------------------------------------------

class PfemSolver:
    def __init__(self, testname, bndno, dt):
        self.testname = testname  # string (name of the module of the fluid model)
        self.bndno = bndno        # phygroup# of the f/s interface
        
        # internal vars
        self.vnods = {}           # dict of interface nodes / prescribed velocities
        self.t1      = 0.0        # last reference time        
        self.t2      = 0.0        # last calculated time
                 
        # loads the python module
        #load(self.testname)         # use toolbox.utilities
        exec("import %s" % self.testname)
        exec("module = %s" % self.testname) # link to Pfem object

        # create an instance of Pfem
        self.pfem = module.getPfem()
        
        # retrieve the f/s boundary and the related nodes
        gr = self.pfem.w.Group(self.pfem.msh, bndno)
        
        # builds a list (dict) of interface nodes
        for e in gr.tag.elems:
            for n in e.nodes:
                no = n.no
                self.vnods[no] = (n)
        
        # Pfem scheme initialization
        self.V = self.pfem.w.DoubleVector()
        self.V0 = self.pfem.w.DoubleVector()
        self.u = self.pfem.w.DoubleVector()
        self.v = self.pfem.w.DoubleVector()
        self.p = self.pfem.w.DoubleVector()
        self.velocity = self.pfem.w.DoubleVector()
        
        self.pfem.scheme.dt = dt
        self.pfem.scheme.init(self.V,self.V0,self.u,self.v,self.p,self.velocity)
        
        self.runOK = True
        
    def run(self, t1, t2):
        """
        calculates one increment from t1 to t2.
        """
        self.pfem.scheme.setNextStep()
        self.runOK = self.pfem.scheme.runOneTimeStep(self.V,self.V0,self.u,self.v,self.p,self.velocity)
        
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
    
    def applyDefo(self, vx, vy, vz, t2):
        """
        prescribes given nodal positions and velocities coming from solid solver to node #no
        """
        if self.pfem.scheme.t < t2:
            self.pfem.scheme.resetNodalPositions()
        
        i = 0
        for no in self.vnods.iterkeys():
            node = self.vnods[no]                 
            node.imposedU = vx[i]
            node.imposedV = vy[i]
            i+=1
            
    def getLoad(self):
        """
        returns the updated loads of the interface nodes.
        """
        fx = np.zeros(len(self.vnods))
        fy = np.zeros(len(self.vnods))
        fz = np.zeros(len(self.vnods))
        
        i = 0
        for no in self.vnods.iterkeys():
            node = self.vnods[no]                 
            fx[i] = -node.Fint.x[0]
            fy[i] = -node.Fint.x[1]
            fz[i] = 0.
            i+=1
        
        return (fx, fy, fz)
    
    def update(self, dt):
        self.pfem.scheme.t+=dt
        self.pfem.scheme.nt+=1
        self.pfem.scheme.updateSolutionVectors(self.V,self.u,self.v,self.p,self.velocity)
        
    def save(self, nt):
        if nt%self.pfem.scheme.savefreq==0:
            self.pfem.scheme.archive()
        self.pfem.scheme.vizu(self.u,self.v,self.p)
        
    def remeshing(self):
        self.pfem.scheme.remeshing(self.V,self.V0,self.p)
        self.pfem.scheme.updateData(self.V0,self.V)
    
    def exit(self):
        """
        Exits the Pfem solver.
        """
        
        self.pfem.gui.save2vtk()
        self.pfem.scheme.exit()
        
        print("***************************** Exit Pfem *****************************")