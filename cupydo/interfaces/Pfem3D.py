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
Python interface between the wrapper of Metafor and CUPyDO.
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
    def __init__(self,p):

        titlePrint('Initializing PFEM3D')
        self.problem = w.getProblem(p['cfdFile'])
        self.interpType = p['interpType']

        # Incompressible or weakly compressible solver

        if 'WC' in self.problem.getID():
            
            self.implicit = False
            self.run = self.runExplicit
            self.maxDivision = 2000

        else:
            
            self.implicit = True
            self.run = self.runImplicit
            self.maxDivision = 10

        # Stores the important objects and variables

        self.FSI = w.VectorInt()
        self.mesh = self.problem.getMesh()
        self.mesh.getNodesIndex('FSInterface',self.FSI)
        self.solver = self.problem.getSolver()
        self.nPhysicalNodes = self.FSI.size()
        self.nNodes = self.FSI.size()

        # Initialize the boundary conditions

        self.BC = list()
        self.dim = self.mesh.getDim()

        for i in self.FSI:

            vector = w.VectorDouble(3)
            self.mesh.getNode(i).setExtState(vector)
            self.BC.append(vector)

        # Save mesh after initializing the BC pointer

        self.prevSolution = w.SolutionData()
        self.problem.copySolution(self.prevSolution)
        self.problem.displayParams()
        self.problem.dump()

        # Store temporary simulation variables

        self.disp = np.zeros((self.nPhysicalNodes,3))
        self.initPos = self.getPosition()
        self.vel = self.getVelocity()
        
        FluidSolver.__init__(self,p)

# Run for implicit integration scheme

    def runImplicit(self,t1,t2):

        print('\nt = {:.5e} - dt = {:.5e}'.format(t2,t2-t1))
        self.problem.loadSolution(self.prevSolution)
        dt = float(t2-t1)
        count = int(1)

        # Main solving loop for the fluid simulation

        while count > 0:
            
            self.solver.setTimeStep(dt)
            if not self.solver.solveOneTimeStep():
                
                dt = float(dt/2)
                count = np.multiply(2,count)
                if dt < (t2-t1)/self.maxDivision: return False
                continue

            count = count-1
        self.__setCurrentState()
        return True

# Run for explicit integration scheme

    def runExplicit(self,t1,t2):

        print('\nt = {:.5e} - dt = {:.5e}'.format(t2,t2-t1))
        self.problem.loadSolution(self.prevSolution)
        iteration = 0

        # Estimate the time step for stability

        self.solver.computeNextDT()
        division = int((t2-t1)/self.solver.getTimeStep())
        if division > self.maxDivision: return False
        dt = (t2-t1)/division

        # Main solving loop for the fluid simulation

        while iteration < division:
    
            iteration += 1
            self.solver.setTimeStep(dt)
            self.solver.solveOneTimeStep()

        self.__setCurrentState()
        return True

# Apply Mechanical Boundary Conditions

    def applyNodalDisplacements(self,dx,dy,dz,haloNodesDisplacements,dt):

        BC = (np.transpose([dx,dy,dz])-self.disp)/dt
        if not self.implicit: BC = 2*(BC-self.vel)/dt

        for i,vector in enumerate(BC):
            for j,val in enumerate(vector): self.BC[i][j] = val

# Apply Thermal Boundary Conditions

    def applyNodalTemperatures(self,Temperature,dt):

        for i,result in enumerate(Temperature):
            self.BC[i][self.dim] = result[0]

# Return Nodal Values

    def getPosition(self):

        result = np.zeros((self.nPhysicalNodes,3))

        for i in range(self.dim):
            for j,k in enumerate(self.FSI):
                result[j,i] = self.mesh.getNode(k).getCoordinate(i)

        return result

    # Computes the nodal velocity vector

    def getVelocity(self):

        result = np.zeros((self.nPhysicalNodes,3))
        
        for i in range(self.dim):
            for j,k in enumerate(self.FSI):
                result[j,i] = self.mesh.getNode(k).getState(i)

        return result

    # Computes the reaction nodal loads

    def __setCurrentState(self):

        result = w.VectorVectorDouble()

        if self.interpType == 'conservative':

            self.solver.computeLoads('FSInterface',self.FSI,result)
            for i in range(self.nNodes):

                self.nodalLoad_X[i] = -result[i][0]
                self.nodalLoad_Y[i] = -result[i][1]
                if self.dim == 3: self.nodalLoad_Z[i] = -result[i][2]

        elif self.dim == 3:

            self.solver.computeStress('FSInterface',self.FSI,result)
            for i in range(self.nNodes):

                self.nodalLoad_XX[i] = result[i][0]
                self.nodalLoad_YY[i] = result[i][1]
                self.nodalLoad_ZZ[i] = result[i][2]
                self.nodalLoad_XY[i] = result[i][3]
                self.nodalLoad_XZ[i] = result[i][4]
                self.nodalLoad_YZ[i] = result[i][5]

        elif self.mesh.isAxiSym():

            self.solver.computeStress('FSInterface',self.FSI,result)
            for i in range(self.nNodes):

                self.nodalLoad_XX[i] = result[i][0]
                self.nodalLoad_YY[i] = result[i][1]
                self.nodalLoad_ZZ[i] = result[i][2]
                self.nodalLoad_XY[i] = result[i][3]

        else:

            self.solver.computeStress('FSInterface',self.FSI,result)
            for i in range(self.nNodes):

                self.nodalLoad_XX[i] = result[i][0]
                self.nodalLoad_YY[i] = result[i][1]
                self.nodalLoad_XY[i] = result[i][2]

# Other Functions

    def update(self,dt):

        self.mesh.remesh(False)
        if self.implicit: self.solver.precomputeMatrix()
        self.problem.copySolution(self.prevSolution)
        self.disp = self.getPosition()-self.initPos
        self.vel = self.getVelocity()

    # Other utilitary functions

    def getNodalIndex(self,index):
        return index

    def getNodalInitialPositions(self):
        return np.transpose(self.initPos)

# Print Functions

    def exit(self):

        self.problem.displayTimeStats()
        titlePrint('Exit PFEM3D')

    # Save te results into a file

    def save(self,_):
        self.problem.dump()
