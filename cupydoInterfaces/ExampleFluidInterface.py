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
# import fluidWrapper
# import math

# Those are mandatory
import numpy as np
from cupydo.genericSolvers import FluidSolver

# ----------------------------------------------------------------------
#  ExampSolver class
# ----------------------------------------------------------------------
               
class ExampSolver(FluidSolver):
    def __init__(self, config): # You are free to add any arguments here
        """
        Des.
        """
        
        print '\n***************************** Initializing Example *****************************'
        

        #self.nNodes =                              # number of nodes (physical + ghost) at the f/s boundary
        #self.nHaloNode =                           # number of ghost nodes at the f/s boundary
        #self.nPhysicalNodes =                      # number of physical nodes at the f/s boundary
        
        FluidSolver.__init__(self)
        self.coreSolver = fluidWrapper.solverDriver(config)
        
        self.initRealTimeData()

    def run(self, t1, t2):
        """
        Des.
        """

        #Run the solver for one iteration, e.g. :
        #self.coreSolver.run()

        self.__setCurrentState()       # use to fill the arrays with nodal values after each run

    def __setCurrentState(self):
        """
        Des.
        """

        #This is an example, you are free to do it your own way
        #for iVertex in range(self.nPhysicalNodes):
            #self.nodalLoad_X[iVertex] = self.coreSolver.getLoadX(iVertex)
            #self.nodalLoad_Y[iVertex] = self.coreSolver.getLoadY(iVertex)
            #self.nodalLoad_Z[iVertex] = self.coreSolver.getLoadZ(iVertex)

    def getNodalInitialPositions(self):
        """
        Des.
        """

        #nodalInitialPos_X = np.zeros(self.nPhysicalNodes)
        #nodalInitialPos_Y = np.zeros(self.nPhysicalNodes)
        #nodalInitialPos_Z = np.zeros(self.nPhysicalNodes)

        #for ii in range(self.nPhysicalNodes):
            #nodalInitialPos_X[ii] =   
            #nodalInitialPos_Y[ii] = 
            #nodalInitialPos_Z[ii] = 

        #return (nodalInitialPos_X, nodalInitialPos_Y, nodalInitialPos_Z)

    def getNodalIndex(self, iVertex):
        """
        Des.
        """

        #no = 

        #return no

    def applyNodalDisplacements(self, dx, dy, dz, dx_nM1, dy_nM1, dz_nM1, haloNodesDisplacements,time):
        """
        Des.
        """

        #This is just an example again
        #for iVertex in range(self.nPhysicalNodes):
            #self.coreSolver.applyNodalDispX(dx[iVertex], iVertex)
            #self.coreSolver.applyNodalDispY(dy[iVertex], iVertex)
            #self.coreSolver.applyNodalDispZ(dz[iVertex], iVertex)

    def update(self, dt):
        """
        Des.
        """

        FluidSolver.update(self)

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

        return

    def initRealTimeData(self):
        """
        Des.
        """
        
        solFile = open('ExampleSolution.ascii', "w")
        solFile.write("Time\tnIter\tValue\n")
        solFile.close()

    def saveRealTimeData(self, time, nFSIIter):
        """
        Des.
        """
        
        solFile = open('ExampleSolution.ascii', "a")
        solFile.write(str(time) + '\t' + str(nFSIIter) + str(1.0) + '\n')
        solFile.close()

    def printRealTimeData(self, time, nFSIIter):
        """
        Des.
        """
        
        toPrint = 'RES-FSI-' + 'ExampleSolution' + ': ' + str(1.0) + '\n'
        print toPrint
    
    def exit(self):
        """
        Des.
        """

        print("***************************** Exit Example solver *****************************")
