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
        self.solverType = self.solver.getType()
        self.mesh = self.problem.getMesh()

        self.thermal = p['thermal']
        self.mechanical = p['mechanical']
        self.interpType = p['interpType']

        # Incompressible or weakly compressible solver

        if self.solverType ==  w.SolverType_Explicit:
            self.max_division = 100

        elif self.solverType == w.SolverType_Implicit:
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

            if self.solverType == w.SolverType_Explicit:
                self.vel = self.__getVelocity()
        
        FluidSolver.__init__(self,p)
        self.initPos = self.__getPosition()
        self.__setCurrentState()

# Calculates One Time Step

    def run(self, t1, t2):

        self.problem.loadSolution(self.prevSolution)
        self.problem.setMinTimeStep((t2-t1)/self.max_division)
        self.problem.setMaxSimTime(t2)

        if self.solverType ==  w.SolverType_Explicit:
            self.solver.computeNextDT()
        
        elif self.solverType == w.SolverType_Implicit:
            self.solver.setTimeStep(t2-t1)

        if self.problem.simulate():

            self.__setCurrentState()
            self.problem.getMesh().savePredictor(True)
            return True

        else: return False

# Apply Mechanical Boundary Conditions

    def applyNodalDisplacements(self, dx, dy, dz, dt, haloNodesDisplacements):

        BC = (np.transpose([dx[0],dy[0],dz[0]])-self.disp)/dt

        if self.solverType == w.SolverType_Explicit:
            BC = 2*(BC-self.vel)/dt

        for i,vector in enumerate(BC):
            for j,val in enumerate(vector): self.BC[i][j] = val

# Apply Thermal Boundary Conditions

    def applyNodalTemperatures(self, Temperature, dt, haloNodesTemperature):

        for i,result in enumerate(Temperature):
            self.BC[i][3] = result

# Return Nodal Values

    def __getPosition(self):
        
        nodes_list = self.problem.getMesh().getNodesList()
        result = np.zeros((self.nPhysicalNodes,3))

        for i in range(self.dim):
            for j,k in enumerate(self.FSI):
                result[j,i] = nodes_list.getCoordinate(k, i)

        return result

# Computes the nodal velocity vector

    def __getVelocity(self):

        nodes_list = self.problem.getMesh().getNodesList()
        result = np.zeros((self.nPhysicalNodes,3))
        
        for i in range(self.dim):
            for j,k in enumerate(self.FSI):
                result[j,i] = nodes_list.getState(k, i)

        return result

# Computes the reaction nodal loads

    def __setCurrentState(self):

        result = w.VectorVectorDouble()

        if self.mechanical: 
            if self.interpType == 'conservative':

                self.solver.computeLoads('FSI', self.FSI,result)
                for i in range(self.nNodes):

                    self.nodalLoad_X[i, 0] = -result[i][0]
                    self.nodalLoad_Y[i, 0] = -result[i][1]
                    if self.dim == 3: self.nodalLoad_Z[i] = -result[i][2]

            elif self.dim == 3:

                self.solver.computeStress('FSI', self.FSI,result)
                for i in range(self.nNodes):

                    self.nodalLoad_XX[i, 0] = result[i][0]
                    self.nodalLoad_YY[i, 0] = result[i][1]
                    self.nodalLoad_ZZ[i, 0] = result[i][2]
                    self.nodalLoad_XY[i, 0] = result[i][3]
                    self.nodalLoad_XZ[i, 0] = result[i][4]
                    self.nodalLoad_YZ[i, 0] = result[i][5]

            elif self.problem.isAxiSymmetric():

                self.solver.computeStress('FSI', self.FSI,result)
                for i in range(self.nNodes):

                    self.nodalLoad_XX[i, 0] = result[i][0]
                    self.nodalLoad_YY[i, 0] = result[i][1]
                    self.nodalLoad_ZZ[i, 0] = result[i][2]
                    self.nodalLoad_XY[i, 0] = result[i][3]

            else:

                self.solver.computeStress('FSI', self.FSI,result)
                for i in range(self.nNodes):

                    self.nodalLoad_XX[i, 0] = result[i][0]
                    self.nodalLoad_YY[i, 0] = result[i][1]
                    self.nodalLoad_XY[i, 0] = result[i][2]

        if self.thermal:
            self.solver.computeHeatFlux('FSI', self.FSI,result)

            for i in range(self.nNodes):

                self.nodalHeatFlux_X[i] = result[i][0]
                self.nodalHeatFlux_Y[i] = result[i][1]
                if self.dim == 3: self.nodalHeatFlux_Z[i] = result[i][2]

# Reset the FSI and Boundary Condition

    def _resetInterfaceBC(self):

        self.BC = list()
        nodes_list = self.problem.getMesh().getNodesList()
        self.mesh.getNodesIndex('FSI', self.FSI)

        for i in self.FSI:

            vector = w.VectorDouble(4)
            nodes_list.setExtState(i, vector)
            self.BC.append(vector)

# Other Functions

    def update(self, dt):

        self.problem.getSolver().updateMesh(True, False)

        # Update the backup and precompute matrices

        self._resetInterfaceBC()
        self.problem.copySolution(self.prevSolution)

        if self.mechanical:
            self.disp = self.__getPosition()-self.initPos

            if self.solverType == w.SolverType_Explicit:
                self.vel = self.__getVelocity()

    def getNodalInitialPositions(self):
        return np.transpose(self.initPos).copy()

# Print Functions
    
    def exit(self):

        self.problem.displayTimeStats()
        titlePrint('Exit PFEM3D')

    def save(self,_):
        self.problem.dump()
