#! /usr/bin/env python3
# -*- coding: utf-8; -*-

'''

Copyright 2018 University of Liège
Copyright 2018 TU Kaiserslautern

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

pyBeamInterface.py
Python interface between the wrapper of pyBeam and CUPyDO.
Authors R. Sanchez - TU Kaiserslautern

'''

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import pyBeamInterface as pyBeam
import math
import numpy as np

# Those are mandatory
import numpy as np
from ..genericSolvers import SolidSolver, SolidAdjointSolver


# ----------------------------------------------------------------------
#  pyBeamSolver class
# ----------------------------------------------------------------------

class pyBeamSolver(SolidSolver):
    """
    pyBeam solver interface.
    """

    def __init__(self, p):
        """
        Initialize the pyBeam solver and all the required interface variables.
        """

        # Initialize pyBeam
        self.initializeSolver(p['csdFile'])

        # Set thickness (Temporary, until we have a config file in place)
        #self.pyBeam.SetThickness(0.02)

        # Initialize solver (will need the config file)
        #self.pyBeam.Initialize()

        # Some options that we should keep just in case
        self.computationType = p['compType']  # computation type : steady (default) or unsteady
        self.nodalLoadsType = p['nodalLoadsType']  # nodal loads type to extract : force (in N, default) or pressure (in Pa)

        # --- Calculate the number of nodes (on each partition) --- #
        self.nNodes = self.pyBeam.nPoint
        self.nHaloNode = 0  # no paralellization required
        self.nPhysicalNodes = self.nNodes  # numbers of nodes at the f/s interface (physical)

        self.nodalInitialPos_X = np.zeros((self.nPhysicalNodes))  # initial position of the f/s interface
        self.nodalInitialPos_Y = np.zeros((self.nPhysicalNodes))
        self.nodalInitialPos_Z = np.zeros((self.nPhysicalNodes))

        # Initialize list to store point indices
        self.pointIndexList = np.zeros(self.nPhysicalNodes, dtype=int)
        self.haloNodesPositionsInit = {}

        self.initializeVariables()

        # --- Initialize the interface position and the nodal loads --- #
        PhysicalIndex = 0
        for iVertex in range(self.nNodes):
            posX, posY, posZ = self.pyBeam.GetInitialCoordinates(iVertex)

            GlobalIndex = iVertex
            self.pointIndexList[PhysicalIndex] = GlobalIndex
            self.nodalInitialPos_X[PhysicalIndex] = posX
            self.nodalInitialPos_Y[PhysicalIndex] = posY
            self.nodalInitialPos_Z[PhysicalIndex] = posZ
            PhysicalIndex += 1

        self.extractors = p['extractors']
        self.__setCurrentState()
        self.initRealTimeData()

        # print("\n -------------------------- SOLID NODES ------------------------------ \n")
        # print(("There is a total of", self.nNodes, "nodes\n"))
        # for iIndex in range(self.nPhysicalNodes):
        #     print((self.pointIndexList[iIndex], self.nodalInitialPos_X[iIndex], self.nodalInitialPos_Y[iIndex], self.nodalInitialPos_Z[iIndex]))

    def initializeSolver(self, confFile):
        self.pyBeam = pyBeam.pyBeamSolver(confFile)

    def initializeVariables(self):
        """
        Initialize variables required by the solver
        """
        SolidSolver.__init__(self)

    def run(self, t1, t2):
        """
        Run one computation of pyBeam.
        """

        if self.computationType == 'unsteady':
            self.__unsteadyRun(t1, t2)
        else:
            self.__steadyRun()

        self.__setCurrentState()
        return True

    def __unsteadyRun(self, t1, t2):
        """
        Run pyBeam on one time step.
        """

        #self.pyBeam.ResetConvergence()
        self.pyBeam.run()

    def __steadyRun(self):
        """
        Run pyBeam up to a converged steady state.
        """

        self.pyBeam.run()

    def __setCurrentState(self):
        """
        Get the nodal (physical) displacements and velocities from pyBeam structural solver.
        """

        PhysicalIndex = 0
        for iVertex in range(self.nNodes):

            #disp = self.SU2.GetDisplacements(self.fluidInterfaceID, iVertex)
            #vel = self.SU2.GetVelocity(self.fluidInterfaceID, iVertex)
            #vel_n = self.SU2.GetVelocity_n(self.fluidInterfaceID, iVertex)

            self.nodalDisp_X[PhysicalIndex], self.nodalDisp_Y[PhysicalIndex], self.nodalDisp_Z[PhysicalIndex] = self.pyBeam.ExtractDisplacements(iVertex)

            # self.nodalVel_X[PhysicalIndex] = vel[0]
            # self.nodalVel_Y[PhysicalIndex] = vel[1]
            # self.nodalVel_Z[PhysicalIndex] = vel[2]

            PhysicalIndex += 1

    def getNodalInitialPositions(self):


        return (self.nodalInitialPos_X, self.nodalInitialPos_Y, self.nodalInitialPos_Z)


    def getNodalDisplacements(self):


        return (self.nodalDisp_X, self.nodalDisp_Y, self.nodalDisp_Z)

    def applyNodalForce(self, load_X, load_Y, load_Z, dt, haloNodesLoads):


        # --- Initialize the interface position and the nodal loads --- #
        PhysicalIndex = 0
        for iVertex in range(self.nNodes):

            LoadX = load_X[PhysicalIndex]
            LoadY = load_Y[PhysicalIndex]
            LoadZ = load_Z[PhysicalIndex]
            PhysicalIndex += 1

            self.pyBeam.SetLoads(iVertex, LoadX, LoadY, LoadZ)

    def update(self):
        """
        Pushes back the current state in the past (previous state) before going to the next time step.
        """

        SolidSolver.update(self)

    def save(self):


        return

    def initRealTimeData(self):


        solFile = open('SolidSolution.ascii', "w")
        solFile.write('{:>12s}   {:>12s}'.format('Time', 'Iteration'))
        for gidx in self.extractors:
            solFile.write('   {:>12s}   {:>12s}   {:>12s}'.format('x_'+str(gidx), 'y_'+str(gidx), 'z_'+str(gidx)))
        solFile.write('\n')
        solFile.close()

    def saveRealTimeData(self, time, nFSIIter):


        solFile = open('SolidSolution.ascii', "a")
        solFile.write("{:>12.6f}   {:>12d}".format(time, nFSIIter))
        for gidx in self.extractors:
            solFile.write('   {:>12.10f}   {:>12.10f}   {:>12.10f}'.format(self.nodalDisp_X[gidx], self.nodalDisp_Y[gidx], self.nodalDisp_Z[gidx]))
        solFile.write('\n')
        solFile.close()

    def printRealTimeData(self, time, nFSIIter):


        toPrint = 'RES-FSI-' + 'pyBeam' + ': ' + str(1.0) + '\n'
        print(toPrint)

    def exit(self):


        print("***************************** Exit pyBeam solver *****************************")

# ----------------------------------------------------------------------
#  pyBeamAdjointSolver class
# ----------------------------------------------------------------------

class pyBeamAdjointSolver(pyBeamSolver, SolidAdjointSolver):
    def initializeSolver(self, confFile):
        self.pyBeam = pyBeam.pyBeamADSolver(confFile)
        self.pyBeam.beam.ReadRestart() # Read restart file and apply displacements

    def initializeVariables(self):
        """
        Initialize variables required by the solver
        """
        SolidAdjointSolver.__init__(self)

    def run(self, t1, t2):
        """
        Run one computation of pyBeam.
        """
        self.__steadyRun()
        self.__setCurrentState()
        return True

    def __steadyRun(self):
        self.pyBeam.RecordSolver()
        self.pyBeam.RunAdjoint()

    def __setCurrentState(self):
        PhysicalIndex = 0
        for iVertex in range(self.nNodes):
            self.nodalDisp_X[PhysicalIndex], self.nodalDisp_Y[PhysicalIndex], self.nodalDisp_Z[PhysicalIndex] = self.pyBeam.ExtractDisplacements(iVertex)
            self.nodalAdjLoad_X[PhysicalIndex], self.nodalAdjLoad_Y[PhysicalIndex], self.nodalAdjLoad_Z[PhysicalIndex] = self.pyBeam.GetLoadAdjoint(iVertex)
            PhysicalIndex += 1

    def applyNodalAdjointDisplacement(self, disp_adj_X, disp_adj_Y, disp_adj_Z, haloNodesDisplacements, dt):
        PhysicalIndex = 0
        for iVertex in range(self.nNodes):
            self.pyBeam.SetDisplacementAdjoint(PhysicalIndex, disp_adj_X[PhysicalIndex], disp_adj_Y[PhysicalIndex], disp_adj_Z[PhysicalIndex])
            PhysicalIndex += 1
