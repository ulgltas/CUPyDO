#! /usr/bin/env python
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

'''

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

# Import whatever you want here, e.g.:
import pythonVLM.VLM_driver as VLM_driver
# import math

# Those are mandatory
import numpy as np
from cupydo.genericSolvers import FluidSolver

# ----------------------------------------------------------------------
#  VLMSolver class
# ----------------------------------------------------------------------
               
class VLMSolver(FluidSolver):
    def __init__(self, infile): # You are free to add any arguments here
        """
        Des.
        """
        
        print '\n***************************** Initializing VLM *****************************'
        
        self.coreSolver = VLM_driver.VLMDriver(infile)
        self.isRun = False
        self.nNodes =   self.coreSolver.data.wing.nvert+ \
                        self.coreSolver.data.flap.nvert+self.coreSolver.data.aileron.nvert  # number of nodes (physical + ghost) at the f/s boundary
        self.nHaloNode = 0                                                                  # number of ghost nodes at the f/s boundary
        self.nPhysicalNodes = self.nNodes-self.nHaloNode                                    # number of physical nodes at the f/s boundary

        FluidSolver.__init__(self)
        
        self.initRealTimeData()

    def run(self, t1, t2):
        """
        Des.
        """

        #Run the solver for one iteration, e.g. :
        if self.isRun:
            self.update(t2-t1)
        else:
            self.isRun = True
        self.coreSolver.run()

        self.__setCurrentState()       # use to fill the arrays with nodal values after each run

    def __setCurrentState(self):
        """
        Des.
        """
        force = np.zeros([self.nPhysicalNodes,3,4])
        weights = [0.375, 0.375, 0.125, 0.125]
        for i in range(self.coreSolver.data.wing.nface):
            indices = self.coreSolver.getVertices(i)
            j = -1
            for index in indices:
                j+=1
                if index>=0:
                    force[index, :, j] = self.coreSolver.getForce(i, weights[j])
        vertexForce = np.sum(force, axis=2)
        #This is an example, you are free to do it your own way
        for iVertex in range(self.nPhysicalNodes):
            self.nodalLoad_X[iVertex] = vertexForce[iVertex, 0]
            self.nodalLoad_Y[iVertex] = vertexForce[iVertex, 1]
            self.nodalLoad_Z[iVertex] = vertexForce[iVertex, 2]

    def getNodalInitialPositions(self):
        """
        Des.
        """

        nodalInitialPos_X = np.zeros(self.nPhysicalNodes)
        nodalInitialPos_Y = np.zeros(self.nPhysicalNodes)
        nodalInitialPos_Z = np.zeros(self.nPhysicalNodes)

        for ii in range(self.nPhysicalNodes):
            nodalInitialPos_X[ii] = self.coreSolver.getX(ii)
            nodalInitialPos_Y[ii] = self.coreSolver.getY(ii)
            nodalInitialPos_Z[ii] = self.coreSolver.getZ(ii)

        return (nodalInitialPos_X, nodalInitialPos_Y, nodalInitialPos_Z)

    def getNodalIndex(self, iVertex):
        """
        Des.
        """
        if iVertex>self.coreSolver.data.wing.nvert+self.coreSolver.data.flap.nvert:
            index = iVertex+20000-self.coreSolver.data.wing.nvert-self.coreSolver.data.flap.nvert
        elif iVertex>self.coreSolver.data.wing.nvert:
            index = iVertex+10000-self.coreSolver.data.wing.nvert
        else:
            index = iVertex
        return index

    def applyNodalDisplacements(self, dx, dy, dz, dx_nM1, dy_nM1, dz_nM1, haloNodesDisplacements,time):
        """
        Des.
        """

        for ii in range(self.coreSolver.m): # For each row of panels
            kk = ii*self.coreSolver.n
            # Current vertex
            self.coreSolver.dX(kk, dx[kk])
            self.coreSolver.dY(kk, dy[kk])
            self.coreSolver.dZ(kk, dz[kk])
            # Current collocation point
            self.coreSolver.dXv(kk, 0.75*dx[kk])
            self.coreSolver.dYv(kk, 0.75*dy[kk])
            self.coreSolver.dZv(kk, 0.75*dz[kk])
            for jj in range(1,self.coreSolver.n-1): # For each panel in a row (except first and last ones)
                kk = ii*self.coreSolver.n+jj
                # Current vertex
                self.coreSolver.dX(kk, dx[kk])
                self.coreSolver.dY(kk, dy[kk])
                self.coreSolver.dZ(kk, dz[kk])
                # Current collocation point
                self.coreSolver.dXv(kk, 0.75*dx[kk])
                self.coreSolver.dYv(kk, 0.75*dy[kk])
                self.coreSolver.dZv(kk, 0.75*dz[kk])
                # Previous collocation point
                self.coreSolver.dXv(kk-1, 0.25*dx[kk])
                self.coreSolver.dYv(kk-1, 0.25*dy[kk])
                self.coreSolver.dZv(kk-1, 0.25*dz[kk])
            kk+=1
            # Current vertex
            self.coreSolver.dX(kk, dx[kk])
            self.coreSolver.dY(kk, dy[kk])
            self.coreSolver.dZ(kk, dz[kk])
            # Previous collocation point
            self.coreSolver.dXv(kk-1, 0.25*dx[kk])
            self.coreSolver.dYv(kk-1, 0.25*dy[kk])
            self.coreSolver.dZv(kk-1, 0.25*dz[kk])
            # Current collocation point
            self.coreSolver.dXv(kk, 1.2*dx[kk])
            self.coreSolver.dYv(kk, 1.2*dy[kk])
            self.coreSolver.dZv(kk, 1.2*dz[kk])
            # Current collocation point, previous vortex
            self.coreSolver.dXv(kk, -0.2*dx[kk-1])
            self.coreSolver.dYv(kk, -0.2*dy[kk-1])
            self.coreSolver.dZv(kk, -0.2*dz[kk-1])

    def update(self, dt):
        """
        Des.
        """

        FluidSolver.update(self, dt)
        self.coreSolver.update()
        #overload here

    def bgsUpdate(self):
        """
        Des.
        """

        #overload here

        return

    def save(self, nt):
        """
        Des.
        """

        #overload here
        self.coreSolver.save()
        return

    def initRealTimeData(self):
        """
        Des.
        """
        
        solFile = open('VLMSolution.ascii', "w")
        solFile.write("Time\tnIter\tValue\n")
        solFile.close()

    def saveRealTimeData(self, time, nFSIIter):
        """
        Des.
        """
        
        solFile = open('VLMSolution.ascii', "a")
        solFile.write(str(time) + '\t' + str(nFSIIter) + str(1.0) + '\n')
        solFile.close()

    def printRealTimeData(self, time, nFSIIter):
        """
        Des.
        """
        
        toPrint = 'RES-FSI-' + 'VLMSolution' + ': ' + str(1.0) + '\n'
        print toPrint
    
    def exit(self):
        """
        Des.
        """

        print("***************************** Exit VLM solver *****************************")
