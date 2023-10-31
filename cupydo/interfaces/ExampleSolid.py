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

'''

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

# Import whatever you want here, e.g.:
# import solidWrapper

# Those are mandatory
import numpy as np
from ..genericSolvers import SolidSolver


# ----------------------------------------------------------------------
#  ExampSolver class
# ----------------------------------------------------------------------
               
class ExampSolver(SolidSolver):
    def __init__(self, config): # You are free to add any arguments here

        
        print('\n***************************** Initializing Example *****************************')
        

        #self.nNodes =                              # number of nodes (physical + ghost) at the f/s boundary
        #self.nHaloNode =                           # number of ghost nodes at the f/s boundary
        #self.nPhysicalNodes =                      # number of physical nodes at the f/s boundary
        
        SolidSolver.__init__(self)
        self.coreSolver = solidWrapper.solverDriver(config)
        
        self.__setCurrentState()           # use to fill the arrays with initial nodal values
        self.initRealTimeData()
        
    def run(self, t1, t2):


        #Run the solver for one iteration, e.g. :
        #self.coreSolver.run()

        self.__setCurrentState()       # use to fill the arrays with nodal values after each run
        return True

    def __setCurrentState(self):


        #This is an example, you are free to do it your own way
        #for ii in range(self.nPhysicalNodes):                   
            #self.nodalDisp_X[ii] = self.coreSolver.getDisplacementX(ii)
            #self.nodalDisp_Y[ii] = self.coreSolver.getDisplacementY(ii)
            #self.nodalDisp_Z[ii] = self.coreSolver.getDisplacementZ(ii)
            #self.nodalVel_X[ii] = self.coreSolver.getVelocityX(ii)  
            #self.nodalVel_Y[ii] = self.coreSolver.getVelocityY(ii)  
            #self.nodalVel_Z[ii] = self.coreSolver.getVelocityZ(ii)  
        return

    def getNodalInitialPositions(self):


        #nodalInitialPos_X = np.zeros(self.nPhysicalNodes)
        #nodalInitialPos_Y = np.zeros(self.nPhysicalNodes)
        #nodalInitialPos_Z = np.zeros(self.nPhysicalNodes)

        #for ii in range(self.nPhysicalNodes):
            #nodalInitialPos_X[ii] =   
            #nodalInitialPos_Y[ii] = 
            #nodalInitialPos_Z[ii] = 

        #return (nodalInitialPos_X, nodalInitialPos_Y, nodalInitialPos_Z)
        return
    

    def applyNodalForce(self, load_X, load_Y, load_Z, dt, haloNodesLoads):


        #This is just an example again
        #for ii in range(self.nPhysicalNodes):
            #self.coreSolver.applyNodalLoadX(load_X[ii], ii)
            #self.coreSolver.applyNodalLoadY(load_Y[ii], ii)
            #self.coreSolver.applyNodalLoadZ(load_Z[ii], ii)
        return

    def applyNodalTemperatures(self, Temperature, dt, haloNodesTemperature):
        return

    def update(self):
        """
        Pushes back the current state in the past (previous state) before going to the next time step.
        """

        SolidSolver.update(self)

        #overload here

    def save(self):


        #overload here

        return

    def initRealTimeData(self):

        
        solFile = open('ExampleSolution.ascii', "w")
        solFile.write("Time\tnIter\tValue\n")
        solFile.close()

    def saveRealTimeData(self, time, nFSIIter):

        
        solFile = open('ExampleSolution.ascii', "a")
        solFile.write(str(time) + '\t' + str(nFSIIter) + str(1.0) + '\n')
        solFile.close()

    def printRealTimeData(self, time, nFSIIter):

        
        toPrint = 'RES-FSI-' + 'ExampleSolution' + ': ' + str(1.0) + '\n'
        print(toPrint)
    
    def exit(self):


        print("***************************** Exit Example solver *****************************")
