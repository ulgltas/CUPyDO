#! /usr/bin/env python
# -*- coding: latin-1; -*-

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

SU2SolidInterface.py
Python interface between the wrapper of SU2 for solid mechanics and CUPyDO.
Authors R. Sanchez - TU Kaiserslautern

'''

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import pysu2ad as pysu2
import math
import numpy as np

# Those are mandatory
import numpy as np
from cupydo.genericSolvers import SolidSolver


# ----------------------------------------------------------------------
#  ExampSolver class
# ----------------------------------------------------------------------

class SU2SolidSolver(SolidSolver):
    """
    SU2 solver interface.
    """

    def __init__(self, confFile, bndno, nDim, computationType, nodalLoadsType, have_MPI, MPIComm=None):
        """
        Initialize the SU2 solver and all the required interface variables.
        """

        # --- Instantiate the structural driver of SU2 --- #
        try:
            print("CGeneralDriver")
            pysu2.CGeneralDriver(confFile, 1, nDim, False, MPIComm)
            print("Goes through")
        except TypeError as exception:
            print('A TypeError occured in pysu2.CSingleZoneDriver : ', exception)
            if have_MPI == True:
                print(
                    'ERROR : You are trying to initialize MPI with a serial build of the wrapper. Please, remove the --parallel option that is incompatible with a serial build.')
            else:
                print(
                    'ERROR : You are trying to launch a computation without initializing MPI but the wrapper has been built in parallel. Please add the --parallel option in order to initialize MPI for the wrapper.')

        allMovingMarkersTags = self.SU2.GetAllMovingMarkersTag()  # list containing the tags of all moving markers
        allMarkersID = self.SU2.GetAllBoundaryMarkers()  # dic : allMarkersID['marker_tag'] = marker_ID
        self.fluidInterfaceID = None  # identification of the f/s boundary, currently limited to one boundary, by default the first tag in allMovingMarkersTags
        if not allMovingMarkersTags:
            # raise Exception('No interface for FSI was defined.')
            self.fluidInterfaceID = None
        elif allMovingMarkersTags:
            if allMovingMarkersTags[0] in allMarkersID.keys():
                self.fluidInterfaceID = allMarkersID[allMovingMarkersTags[0]]
            else:
                raise Exception("Moving and CHT markes have to be the same.")

        self.computationType = computationType  # computation type : steady (default) or unsteady
        self.nodalLoadsType = nodalLoadsType  # nodal loads type to extract : force (in N, default) or pressure (in Pa)

        # --- Calculate the number of nodes (on each partition) --- #
        self.nNodes = 0
        self.nHaloNode = 0
        self.nPhysicalNodes = 0
        if self.fluidInterfaceID != None:
            self.nNodes = self.SU2.GetNumberVertices(
                self.fluidInterfaceID)  # numbers of nodes at the f/s interface (halo+physical)
            self.nHaloNode = self.SU2.GetNumberHaloVertices(
                self.fluidInterfaceID)  # numbers of nodes at the f/s interface (halo)
        self.nPhysicalNodes = self.nNodes - self.nHaloNode  # numbers of nodes at the f/s interface (physical)

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
            posX = self.SU2.GetVertexCoordX(self.fluidInterfaceID, iVertex)
            posY = self.SU2.GetVertexCoordY(self.fluidInterfaceID, iVertex)
            posZ = self.SU2.GetVertexCoordZ(self.fluidInterfaceID, iVertex)
            if self.SU2.IsAHaloNode(self.fluidInterfaceID, iVertex):
                GlobalIndex = self.SU2.GetVertexGlobalIndex(self.fluidInterfaceID, iVertex)
                self.haloNodeList[GlobalIndex] = iVertex
                self.haloNodesPositionsInit[GlobalIndex] = (posX, posY, posZ)
            else:
                GlobalIndex = self.SU2.GetVertexGlobalIndex(self.fluidInterfaceID, iVertex)
                self.pointIndexList[PhysicalIndex] = GlobalIndex
                self.SU2.ComputeVertexForces(self.fluidInterfaceID, iVertex)
                Fx = self.SU2.GetVertexForceX(self.fluidInterfaceID, iVertex)
                Fy = self.SU2.GetVertexForceY(self.fluidInterfaceID, iVertex)
                Fz = self.SU2.GetVertexForceZ(self.fluidInterfaceID, iVertex)
                Temp = self.SU2.GetVertexTemperature(self.fluidInterfaceID, iVertex)
                self.nodalInitialPos_X[PhysicalIndex] = posX
                self.nodalInitialPos_Y[PhysicalIndex] = posY
                self.nodalInitialPos_Z[PhysicalIndex] = posZ
                self.nodalTemperature[PhysicalIndex] = Temp
                PhysicalIndex += 1

        self.initRealTimeData()

        print("\n -------------------------- SOLID NODES ------------------------------ \n")
        print("There is a total of", self.nNodes, "nodes\n")
        for iIndex in range(self.nPhysicalNodes):
            print(self.pointIndexList[iIndex], self.nodalInitialPos_X[iIndex], self.nodalInitialPos_Y[iIndex], self.nodalInitialPos_Z[iIndex])

    def run(self, t1, t2):
        """
        Run one computation of SU2.
        """

        if self.computationType == 'unsteady':
            self.__unsteadyRun(t1, t2)
        else:
            self.__steadyRun()

        self.__setCurrentState()

    def __unsteadyRun(self, t1, t2):
        """
        Run SU2 on one time step.
        """

        self.SU2.ResetConvergence()
        self.SU2.Run()

    def __steadyRun(self):
        """
        Run SU2 up to a converged steady state.
        """

        #self.SU2.ResetConvergence()
        self.SU2.Run()

    def __setCurrentState(self):
        """
        Get the nodal (physical) displacements and velocities from SU2 structural solver.
        """

        PhysicalIndex = 0
        for iVertex in range(self.nNodes):
            # identify the halo nodes and ignore their nodal loads
            halo = self.SU2.ComputeVertexForces(self.fluidInterfaceID, iVertex)
            self.SU2.ComputeVertexHeatFluxes(self.fluidInterfaceID, iVertex)
            if halo == False:

                disp = self.SU2.GetDisplacements(self.fluidInterfaceID, iVertex)
                vel = self.SU2.GetVelocity(self.fluidInterfaceID, iVertex)
                vel_n = self.SU2.GetVelocity_n(self.fluidInterfaceID, iVertex)

                self.nodalDisp_X[PhysicalIndex] = disp[0]
                self.nodalDisp_Y[PhysicalIndex] = disp[1]
                self.nodalDisp_Z[PhysicalIndex] = disp[2]

                self.nodalVel_X[PhysicalIndex] = vel[0]
                self.nodalVel_Y[PhysicalIndex] = vel[1]
                self.nodalVel_Z[PhysicalIndex] = vel[2]

                self.nodalVel_XNm1[PhysicalIndex] = vel_n[0]
                self.nodalVel_YNm1[PhysicalIndex] = vel_n[1]
                self.nodalVel_ZNm1[PhysicalIndex] = vel_n[2]

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

    def applyNodalLoads(self, load_X, load_Y, load_Z, val_time):
        """
        Des.
        """

        # --- Initialize the interface position and the nodal loads --- #
        PhysicalIndex = 0
        for iVertex in range(self.nNodes):
            GlobalIndex = self.SU2.GetVertexGlobalIndex(self.fluidInterfaceID, iVertex)
            # in case of halo node, use haloNodesDisplacements with global fluid indexing
            if GlobalIndex in self.haloNodeList.keys():
                # Temporarily support only single core
                LoadX = 0.0
                LoadY = 0.0
                LoadZ = 0.0
            else:
                LoadX = load_X[PhysicalIndex]
                LoadY = load_Y[PhysicalIndex]
                LoadZ = load_Z[PhysicalIndex]
                PhysicalIndex += 1
            self.SU2.SetLoads(self.fluidInterfaceID, iVertex, GlobalIndex, LoadX, LoadY, LoadZ)

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
        print
        toPrint

    def exit(self):
        """
        Des.
        """
        self.SU2.Output(2) # Temporary hack
        self.SU2.Postprocessing()
        print("***************************** Exit Example solver *****************************")
