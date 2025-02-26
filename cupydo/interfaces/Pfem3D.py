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

Pfem3D.py
Python interface between the wrapper of PFEM3D and CUPyDO.
Authors M. Lacroix, S. Février

'''

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

from ..genericSolvers import FluidSolver
from ..utilities import titlePrint
import pfem3Dw as w
import numpy as np

# ----------------------------------------------------------------------
#  Pfem3D solver interface class
# ----------------------------------------------------------------------

class Pfem3D(FluidSolver):
    def __init__(self, p):

        titlePrint('Initializing PFEM3D')
        
        self.problem = w.getProblem(p['cfdFile'])
        self.solver = self.problem.getSolver()
        self.mesh = self.problem.getMesh()

        self.thermal = p['thermal']
        self.mechanical = p['mechanical']
        self.interpType = p['interpType']

        # Incompressible or weakly compressible solver

        if 'WC' in self.problem.getID():

            self.WC = True
            self.max_division = 100

        else:

            self.WC = False
            self.max_division = 1

        # Initialize the boundary conditions

        self.BC = list()
        self.FSI = w.VectorInt()
        self._resetInterfaceBC()

        # Stores the dimension and number of interface nodes

        self.dim = self.mesh.getDim()
        self.nNodes = self.FSI.size()
        self.nPhysicalNodes = self.FSI.size()

        # Save mesh after initializing the BC pointer

        self.prevSolution = w.SolutionData()
        self.problem.copySolution(self.prevSolution)
        self.problem.displayParams()
        self.problem.dump()

        # Store temporary simulation variables

        if self.mechanical:
            self.disp = np.zeros((self.nPhysicalNodes,3))
            self.vel = self.__getVelocity()
        
        FluidSolver.__init__(self,p)
        self.initPos = self.__getPosition()
        self.__setCurrentState()

# Calculates One Time Step

    def run(self, t1, t2):

        self.problem.loadSolution(self.prevSolution)
        self.problem.setMinTimeStep((t2-t1)/self.max_division)
        self.problem.setMaxSimTime(t2)

        if self.WC: self.solver.computeNextDT()
        else: self.solver.setTimeStep(t2-t1)

        if self.problem.simulate():
            self.__setCurrentState()
            return True

        else: return False

# Apply Mechanical Boundary Conditions

    def applyNodalDisplacements(self, dx, dy, dz, dt, haloNodesDisplacements):

        BC = (np.transpose([dx,dy,dz])-self.disp)/dt
        if self.WC: BC = 2*(BC-self.vel)/dt

        for i,vector in enumerate(BC):
            for j,val in enumerate(vector): self.BC[i][j] = val

# Apply Thermal Boundary Conditions

    def applyNodalTemperatures(self, Temperature, dt, haloNodesTemperature):

        for i,result in enumerate(Temperature):
            self.BC[i][3] = result

# Return Nodal Values

    def __getPosition(self):

        result = np.zeros((self.nPhysicalNodes,3))

        for i in range(self.dim):
            for j,k in enumerate(self.FSI):
                result[j,i] = self.mesh.getNode(k).getCoordinate(i)

        return result

# Computes the nodal velocity vector

    def __getVelocity(self):

        result = np.zeros((self.nPhysicalNodes,3))
        
        for i in range(self.dim):
            for j,k in enumerate(self.FSI):
                result[j,i] = self.mesh.getNode(k).getState(i)

        return result

# Computes the reaction nodal loads

    def __setCurrentState(self):

        result = w.VectorVectorDouble()

        if self.mechanical: 
            if self.interpType == 'conservative':

                self.solver.computeLoads('FSInterface', self.FSI,result)
                for i in range(self.nNodes):

                    self.nodalLoad_X[i] = -result[i][0]
                    self.nodalLoad_Y[i] = -result[i][1]
                    if self.dim == 3: self.nodalLoad_Z[i] = -result[i][2]

            elif self.dim == 3:

                self.solver.computeStress('FSInterface', self.FSI,result)
                for i in range(self.nNodes):

                    self.nodalLoad_XX[i] = result[i][0]
                    self.nodalLoad_YY[i] = result[i][1]
                    self.nodalLoad_ZZ[i] = result[i][2]
                    self.nodalLoad_XY[i] = result[i][3]
                    self.nodalLoad_XZ[i] = result[i][4]
                    self.nodalLoad_YZ[i] = result[i][5]

            elif self.mesh.isAxiSym():

                self.solver.computeStress('FSInterface', self.FSI,result)
                for i in range(self.nNodes):

                    self.nodalLoad_XX[i] = result[i][0]
                    self.nodalLoad_YY[i] = result[i][1]
                    self.nodalLoad_ZZ[i] = result[i][2]
                    self.nodalLoad_XY[i] = result[i][3]

            else:

                self.solver.computeStress('FSInterface', self.FSI,result)
                for i in range(self.nNodes):

                    self.nodalLoad_XX[i] = result[i][0]
                    self.nodalLoad_YY[i] = result[i][1]
                    self.nodalLoad_XY[i] = result[i][2]

        if self.thermal:
            self.solver.computeHeatFlux('FSInterface', self.FSI,result)

            for i in range(self.nNodes):

                self.nodalHeatFlux_X[i] = result[i][0]
                self.nodalHeatFlux_Y[i] = result[i][1]
                if self.dim == 3: self.nodalHeatFlux_Z[i] = result[i][2]

# Reset the FSI and Boundary Condition

    def _resetInterfaceBC(self):

        self.BC = list()
        self.mesh.getNodesIndex('FSInterface', self.FSI)

        for i in self.FSI:

            vector = w.VectorDouble(4)
            self.mesh.getNode(i).setExtState(vector)
            self.BC.append(vector)

# Other Functions

    def update(self, dt):

        self.mesh.remesh(verboseOutput = False)
        self._resetInterfaceBC()

        # Update the backup and precompute matrices

        if not self.WC: self.solver.precomputeMatrix()
        self.problem.copySolution(self.prevSolution)

        self.disp = self.__getPosition()-self.initPos
        self.vel = self.__getVelocity()

    def getNodalInitialPositions(self):
        return np.transpose(self.initPos).copy()

# Print Functions
    
    def exit(self):

        self.problem.displayTimeStats()
        titlePrint('Exit PFEM3D')

    def save(self,_):
        self.problem.dump()
