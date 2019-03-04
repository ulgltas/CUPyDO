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
    def __init__(self, _case, _nthreads):       
        # load the python module
        case = __import__(_case)
        
        # load the case (contains flow and some of its objects)
        # TODO: check what needs to be passed to getFlow
        self.flow = case.getFlow()

        # get the f/s boundary
        self.boundary = self.flow.boundary

        # number of nodes
        self.nNodes = self.boundary.nodes.size()
        self.nHaloNode = 0
        self.nPhysicalNodes = self.nNodes - self.nHaloNode

        # initial nodal position
        self.nodalInitPosX, self.nodalInitPosY, self.nodalInitPosZ = self.getNodalInitialPositions()

        # nodal load
        self.nodalLoad_X = np.zeros(self.nPhysicalNodes)
        self.nodalLoad_Y = np.zeros(self.nPhysicalNodes)
        self.nodalLoad_Z = np.zeros(self.nPhysicalNodes)

        # initialize
        self.flow.solver.nthreads = _nthreads
        self.flow.mshDef.nthreads = _nthreads
        self.exeOK = True
        FluidSolver.__init__(self)
        
    def run(self, t1, t2):
        """
        Run the solver for one steady (time) iteration.
        """

        self.exeOK = self.flow.solver.run()
        self.__setCurrentState()
    
    def __setCurrentState(self):
        """
        Compute nodal forces from nodal Pressure coefficient
        """

        # integrate Cp at element
        cpiE = self.boundary.integrate(self.flow.solver.phi, self.flow.fCp)
        # transfer integrated Cp from elements to nodes
        cfN = self.boundary.transfer(cpiE)
        i = 0
        for n in self.boundary.nodes:
            self.nodalLoad_X[i] = -self.flow.dynP * cfN[i][0]
            self.nodalLoad_Y[i] = -self.flow.dynP * cfN[i][1]
            self.nodalLoad_Z[i] = -self.flow.dynP * cfN[i][2]
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

    def applyNodalDisplacements(self, dx, dy, dz, dx_nM1, dy_nM1, dz_nM1, haloNodesDisplacements, time):
        """
        Apply displacements coming from solid solver to f/s interface after saving
        """

        self.flow.mshDef.savePos()
        for i in range(self.boundary.nodes.size()):
            self.boundary.nodes[i].pos.x[0] = self.nodalInitPosX[i] + dx[i]
            self.boundary.nodes[i].pos.x[1] = self.nodalInitPosY[i] + dy[i]
            self.boundary.nodes[i].pos.x[2] = self.nodalInitPosZ[i] + dz[i]

    def meshUpdate(self, nt):
        """
        Deform the mesh using linear elasticity equations
        """
        self.flow.mshDef.deform()
        
    def update(self, dt):
        """
        TODO
        """
        
    def save(self, nt):
        """
        Save data on disk at each converged timestep
        """
        self.flow.solver.save(nt)
        self.flow.msh.save(self.flow.msh.name + "_" + str(nt) + ".msh")

    def initRealTimeData(self):
        """
        Initialize history file
        """
        histFile = open('FlowHistory.dat', "w")
        histFile.write("{0:>12s}   {1:>12s}   {2:>12s}   {3:>12s}\n".format("Time", "FSI_Iter", "C_Lift", "C_Drag"))
        histFile.close()

    def saveRealTimeData(self, time, nFSIIter):
        """
        Save data in history file at each fsi iteration
        """
        histFile = open('FlowHistory.dat', "a")
        histFile.write("{0:12.6f}   {1:12d}   {2:12.6f}   {3:12.6f}\n".format(time, nFSIIter, self.flow.solver.Cl, self.flow.solver.Cd))
        histFile.close()

    def printRealTimeData(self, time, nFSIIter):
        """
        Print data on screen...
        """
        
        #print toPrint
    
    def exit(self):
        """
        Exit the Flow solver.
        """
