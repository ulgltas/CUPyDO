#! /usr/bin/env python
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
from ..genericSolvers import SolidSolver


# ----------------------------------------------------------------------
#  pyBeamSolver class
# ----------------------------------------------------------------------

class pyBeamSolver(SolidSolver):
    """
    pyBeam solver interface.
    """

    def __init__(self, confFile, nDim, computationType, nodalLoadsType, extractors):
        """
        Initialize the pyBeam solver and all the required interface variables.
        """

        # Initialize pyBeam
        self.pyBeam = pyBeam.pyBeamSolver(confFile)

        # Set thickness (Temporary, until we have a config file in place)
        #self.pyBeam.SetThickness(0.02)

        # Initialize solver (will need the config file)
        #self.pyBeam.Initialize()

        # Some options that we should keep just in case
        self.computationType = computationType  # computation type : steady (default) or unsteady
        self.nodalLoadsType = nodalLoadsType  # nodal loads type to extract : force (in N, default) or pressure (in Pa)

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

        SolidSolver.__init__(self)

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

        self.extractors = extractors

        self.initRealTimeData()

        # print("\n -------------------------- SOLID NODES ------------------------------ \n")
        # print(("There is a total of", self.nNodes, "nodes\n"))
        # for iIndex in range(self.nPhysicalNodes):
        #     print((self.pointIndexList[iIndex], self.nodalInitialPos_X[iIndex], self.nodalInitialPos_Y[iIndex], self.nodalInitialPos_Z[iIndex]))

    def run(self, t1, t2):
        """
        Run one computation of pyBeam.
        """

        if self.computationType == 'unsteady':
            self.__unsteadyRun(t1, t2)
        else:
            self.__steadyRun()

        self.__setCurrentState()

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
            #
            # self.nodalVel_XNm1[PhysicalIndex] = vel_n[0]
            # self.nodalVel_YNm1[PhysicalIndex] = vel_n[1]
            # self.nodalVel_ZNm1[PhysicalIndex] = vel_n[2]

            PhysicalIndex += 1

    def getNodalInitialPositions(self):
        """
        Description.
        """

        return (self.nodalInitialPos_X, self.nodalInitialPos_Y, self.nodalInitialPos_Z)

    def getNodalIndex(self, iVertex):
        """
        Returns the index (identifier) of the iVertex^th interface node.
        """

        # no =

        # return no

    def getNodalDisplacements(self):
        """
        Des.
        """

        return (self.nodalDisp_X, self.nodalDisp_Y, self.nodalDisp_Z)

    def applyNodalLoads(self, load_X, load_Y, load_Z, val_time, haloNodesLoads = {}):
        """
        Des.
        """

        # --- Initialize the interface position and the nodal loads --- #
        PhysicalIndex = 0
        for iVertex in range(self.nNodes):

            LoadX = load_X[PhysicalIndex]
            LoadY = load_Y[PhysicalIndex]
            LoadZ = load_Z[PhysicalIndex]
            print((LoadX, LoadY, LoadZ))
            PhysicalIndex += 1

            self.pyBeam.SetLoads(iVertex, LoadX, LoadY, LoadZ)

    def update(self):
        """
        Pushes back the current state in the past (previous state) before going to the next time step.
        """

        SolidSolver.update(self)

        # overload here

    def bgsUpdate(self):
        """
        Des.
        """

        # overload here

        return

    def save(self):
        """
        Des.
        """

        # overload here

        return

    def initRealTimeData(self):
        """
        Des.
        """

        solFile = open('SolidSolution.ascii', "w")
        solFile.write('{:>12s}   {:>12s}'.format('Time', 'Iteration'))
        for gidx in self.extractors:
            solFile.write('   {:>12s}   {:>12s}   {:>12s}'.format('x_'+str(gidx), 'y_'+str(gidx), 'z_'+str(gidx)))
        solFile.write('\n')
        solFile.close()

    def saveRealTimeData(self, time, nFSIIter):
        """
        Des.
        """

        solFile = open('SolidSolution.ascii', "a")
        solFile.write("{:>12.6f}   {:>12d}".format(time, nFSIIter))
        for gidx in self.extractors:
            solFile.write('   {:>12.10f}   {:>12.10f}   {:>12.10f}'.format(self.nodalDisp_X[gidx], self.nodalDisp_Y[gidx], self.nodalDisp_Z[gidx]))
        solFile.write('\n')
        solFile.close()

    def printRealTimeData(self, time, nFSIIter):
        """
        Des.
        """

        toPrint = 'RES-FSI-' + 'pyBeam' + ': ' + str(1.0) + '\n'
        print(toPrint)

    def exit(self):
        """
        Des.
        """

        print("***************************** Exit pyBeam solver *****************************")
