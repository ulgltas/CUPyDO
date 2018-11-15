#!/usr/bin/env python
# -*- coding: utf-8; -*-

"""
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

@file FlowInterface.py
@brief Python interface between Flow and CUPyDO.
@author A. Crovato
"""

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import os, os.path, sys, time, string

import math
import numpy as np
from cupydo.genericSolvers import FluidSolver

# ----------------------------------------------------------------------
#  FlowSolver class
# ----------------------------------------------------------------------

class Flow(FluidSolver):
    def __init__(self, _case):       
        # load the python module
        import _case as case
        
        # load the case (contains flow and some of its objects)
        self.flow = case.getFlow()

        # get the f/s boundary
        self.boundary = self.flow.boundary

        # number of nodes
        self.nNodes = self.boundary.nodes.size()
        self.nHaloNode = 0
        self.nPhysicalNodes = self.nNodes - self.nHaloNode

        # initialize
        self.exeOK = True
        FluidSolver.__init__(self)
        # TODO should nodal positions and loads be initialized here too?
        
    def run(self, t1, t2):
        """
        Run the solver for one steady (time) iteration.
        """
        self.exeOK = self.flow.solver.execute()
        self.__setCurrentState()
    
    def __setCurrentState(self):
        """
        Compute nodal forces from nodal Pressure coefficient
        """
        
        # elemental forces: integrate pressure from nodal potential
        p = 0.5 # dynamic pressure
        if  M_inf == 0:
            CfE = airfoil.integrate(self.flow.solver.phi, self.flow.flow.Fun0EleCpL())
        else:
            CfE = airfoil.integrate(self.flow.solver.phi, self.flow.flow.Fun0EleCp(self.flow.gamma, self.flow.M_inf))
        fxE = {}
        fyE = {}
        fzE = {}
        i = 0
        for e in self.boundary.tag.elems:
            fxE[e] = p * cfE[i] * e.normal().x[0]
            fyE[e] = p * cfE[i] * e.normal().x[1]
            fzE[e] = p * cfE[i] * e.normal().x[2]
            i += 1

        # nodal forces: average forces from elements
        i = 0
        for n in self.boundary.nodes:
            for e in self.boundary.neMap[n]:
                fx += fxE.get(e)
                fy += fyE.get(e)
                fz += fzE.get(e)
            self.nodalLoad_X[i] = fx / self.boundary.neMap[n].size()
            self.nodalLoad_Y[i] = fy / self.boundary.neMap[n].size()
            self.nodalLoad_Z[i] = fz / self.boundary.neMap[n].size()
            i += 1

    def getNodalInitialPositions(self):
        """
        Get the initial position of each node
        """
        
        x0 = np.zeros(self.nPhysicalNodes)
        y0 = np.zeros(self.nPhysicalNodes)
        z0 = np.zeros(self.nPhysicalNodes)
        for i in range(self.boundary.nodes.size()):
            n = self.boundary.nodes[i]               
            x0[i] = n.pos.x[0]
            y0[i] = n.pos.x[1]
            z0[i] = n.pos.x[2]

        return (x0, y0, z0)

    def getNodalIndex(self, iVertex):
        """
        Get index of each node
        """

        no = self.boundary.nodes[iVertex].no
        return no

    def applyNodalDisplacements(self, disp_X, disp_Y, disp_Z):
        """
        Apply displacements coming from solid solver to f/s interface
        """

        for i in range(self.boundary.nodes.size()):
            self.boundary.nodes[i].pos.x[0] += disp_X[i]
            self.boundary.nodes[i].pos.x[1] += disp_Y[i]
            self.boundary.nodes[i].pos.x[2] += disp_Z[i]

    def remeshing(self):
        """
        TODO Lagrangian remeshing
        """

    def fakeSolidSolver(self):
        """
        Dummy solid solver for testing
        Apply dummy displacement
        """

        dxD = np.zeros(self.nPhysicalNodes)
        dyD = np.zeros(self.nPhysicalNodes)
        dzD = np.zeros(self.nPhysicalNodes)
        self.applyNodalDisplacements(self, dxD, dyD, dzD) # TODO: check 'self' parameter           
        
    def update(self):
        """
        TODO
        """
        
    def save(self):
        """
        TODO
        """
        self.flow.solver.finalize()

    def initRealTimeData(self):
        """
        TODO
        """

    def saveRealTimeData(self, time, nFSIIter):
        """
        TODO
        """

    def printRealTimeData(self, time, nFSIIter):
        """
        TODO
        """
        
        #print toPrint
    
    def exit(self):
        """
        Exit the Flow solver.
        """
        
        print("***************************** Exit Flow *****************************")
