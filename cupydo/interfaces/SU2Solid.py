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

SU2Solid.py
Python interface between the wrapper of SU2 for solid mechanics and CUPyDO.
Authors R. Sanchez - TU Kaiserslautern

'''
from __future__ import print_function
from __future__ import absolute_import

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

from builtins import str
from builtins import range
import pysu2ad as pysu2
import math
import numpy as np

# Those are mandatory
import numpy as np
from ..genericSolvers import SolidSolver, SolidAdjointSolver


# ----------------------------------------------------------------------
#  ExampSolver class
# ----------------------------------------------------------------------

class SU2SolidSolver(SolidSolver):
    """
    SU2 solver interface.
    """

    def __init__(self, confFile, nDim, computationType, nodalLoadsType, have_MPI, MPIComm=None):
        """
        Initialize the SU2 solver and all the required interface variables.
        """

        self.initializeSolver(confFile, have_MPI, MPIComm)
        allFluidLoadMarkersTags = self.SU2.GetAllFluidLoadMarkersTag()  # list containing the tags of all fluid load markers
        allMarkersID = self.SU2.GetAllBoundaryMarkers()  # dic : allMarkersID['marker_tag'] = marker_ID
        self.solidInterfaceID = None  # identification of the f/s boundary, currently limited to one boundary, by default the first tag in allFluidLoadMarkersTags
        if not allFluidLoadMarkersTags:
            # raise Exception('No interface for FSI was defined.')
            self.solidInterfaceID = None
        elif allFluidLoadMarkersTags:
            if allFluidLoadMarkersTags[0] in list(allMarkersID.keys()):
                self.solidInterfaceID = allMarkersID[allFluidLoadMarkersTags[0]]
            else:
                raise Exception("Moving and CHT markes have to be the same.")

        self.computationType = computationType  # computation type : steady (default) or unsteady
        self.nodalLoadsType = nodalLoadsType  # nodal loads type to extract : force (in N, default) or pressure (in Pa)

        # --- Calculate the number of nodes (on each partition) --- #
        self.nNodes = 0
        self.nHaloNode = 0
        self.nPhysicalNodes = 0
        if self.solidInterfaceID != None:
            self.nNodes = self.SU2.GetNumberVertices(
                self.solidInterfaceID)  # numbers of nodes at the f/s interface (halo+physical)
            self.nHaloNode = self.SU2.GetNumberHaloVertices(
                self.solidInterfaceID)  # numbers of nodes at the f/s interface (halo)
        self.nPhysicalNodes = self.nNodes - self.nHaloNode  # numbers of nodes at the f/s interface (physical)

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
            posX = self.SU2.GetVertexCoordX(self.solidInterfaceID, iVertex)
            posY = self.SU2.GetVertexCoordY(self.solidInterfaceID, iVertex)
            posZ = self.SU2.GetVertexCoordZ(self.solidInterfaceID, iVertex)
            if self.SU2.IsAHaloNode(self.solidInterfaceID, iVertex):
                GlobalIndex = self.SU2.GetVertexGlobalIndex(self.solidInterfaceID, iVertex)
                self.haloNodeList[GlobalIndex] = iVertex
                self.haloNodesPositionsInit[GlobalIndex] = (posX, posY, posZ)
            else:
                GlobalIndex = self.SU2.GetVertexGlobalIndex(self.solidInterfaceID, iVertex)
                self.pointIndexList[PhysicalIndex] = GlobalIndex
                self.SU2.ComputeVertexForces(self.solidInterfaceID, iVertex)
                Fx = self.SU2.GetVertexForceX(self.solidInterfaceID, iVertex)
                Fy = self.SU2.GetVertexForceY(self.solidInterfaceID, iVertex)
                Fz = self.SU2.GetVertexForceZ(self.solidInterfaceID, iVertex)
                Temp = self.SU2.GetVertexTemperature(self.solidInterfaceID, iVertex)
                self.nodalInitialPos_X[PhysicalIndex] = posX
                self.nodalInitialPos_Y[PhysicalIndex] = posY
                self.nodalInitialPos_Z[PhysicalIndex] = posZ
                self.nodalTemperature[PhysicalIndex] = Temp
                PhysicalIndex += 1

        self.initRealTimeData()

        # print("\n -------------------------- SOLID NODES ------------------------------ \n")
        # print(("There is a total of", self.nNodes, "nodes\n"))
        # for iIndex in range(self.nPhysicalNodes):
        #     print((self.pointIndexList[iIndex], self.nodalInitialPos_X[iIndex], self.nodalInitialPos_Y[iIndex], self.nodalInitialPos_Z[iIndex]))

    def initializeSolver(self, confFile, have_MPI, MPIComm):
        # --- Instantiate the single zone driver of SU2 --- #
        # @todo [Adrien Crovato ]Change CFluidDriver constructor
        # as of SU2-6.1.0, a new way of handling periodic boundary conditon has been implemented
        # Consequently CDriver(config, nZone, nDim, MPIComm) changed to CFluidDriver(config, nZone, nDim, val_periodic, MPIComm)
        # Since periodic BC are not used yet in CUPyDO, I just adapted the constructor. This will have to be changed...
        try:
            self.SU2 = pysu2.CSinglezoneDriver(confFile, 1, MPIComm)
        except TypeError as exception:
            print(('A TypeError occured in pysu2.CFluidDriver : ',exception))
            if have_MPI == True:
                print('ERROR : You are trying to initialize MPI with a serial build of the wrapper. Please, remove the --parallel option that is incompatible with a serial build.')
            else:
                print('ERROR : You are trying to launch a computation without initializing MPI but the wrapper has been built in parallel. Please add the --parallel option in order to initialize MPI for the wrapper.')
            
    def initializeVariables(self):
        """
        Initialize variables required by the solver
        """
        SolidSolver.__init__(self)

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

        self.SU2.ResetConvergence()
        self.SU2.Preprocess(0)
        self.SU2.Run()
        StopIntegration = self.SU2.Monitor(0)
        self.SU2.Output(0)

    def setInitialDisplacements(self):
        """
        Set initial displacements
        """
        self.__setCurrentState()

    def __setCurrentState(self):
        """
        Get the nodal (physical) displacements and velocities from SU2 structural solver.
        """

        PhysicalIndex = 0
        for iVertex in range(self.nNodes):
            # identify the halo nodes and ignore their nodal loads
            halo = self.SU2.ComputeVertexForces(self.solidInterfaceID, iVertex)
            self.SU2.ComputeVertexHeatFluxes(self.solidInterfaceID, iVertex)
            if halo == False:

                disp = self.SU2.GetFEA_Displacements(self.solidInterfaceID, iVertex)
                vel = self.SU2.GetFEA_Velocity(self.solidInterfaceID, iVertex)
                vel_n = self.SU2.GetFEA_Velocity_n(self.solidInterfaceID, iVertex)

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

        return self.SU2.GetVertexGlobalIndex(self.solidInterfaceID, int(iVertex))

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
            GlobalIndex = self.SU2.GetVertexGlobalIndex(self.solidInterfaceID, iVertex)
            # in case of halo node, use haloNodesDisplacements with global fluid indexing
            if GlobalIndex in list(self.haloNodeList.keys()):
                # Temporarily support only single core
                LoadX = 0.0
                LoadY = 0.0
                LoadZ = 0.0
            else:
                LoadX = load_X[PhysicalIndex]
                LoadY = load_Y[PhysicalIndex]
                LoadZ = load_Z[PhysicalIndex]
                PhysicalIndex += 1
            self.SU2.SetFEA_Loads(self.solidInterfaceID, iVertex, LoadX, LoadY, LoadZ)

    def update(self):
        """
        Pushes back the current state in the past (previous state) before going to the next time step.
        """

        self.SU2.Update()

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

        stopComp = self.SU2.Monitor(1)
        self.SU2.Output(1)

        return stopComp

    def initRealTimeData(self):
        """
        Des.
        """

        solFile = open('SolidSolution.ascii', "w")
        solFile.write("Time\tnIter\tY_LE\tY_TE\n")
        solFile.close()

    def saveRealTimeData(self, time, nFSIIter):
        """
        Des.
        """

        solFile = open('SolidSolution.ascii', "a")
        solFile.write("{}\t{}\t{}\t{}\n".format(time, nFSIIter, self.nodalDisp_Y[-1], self.nodalDisp_Y[2]))
        solFile.close()

    def printRealTimeData(self, time, nFSIIter):
        """
        Des.
        """

        toPrint = 'RES-FSI-' + 'SU2SolidSolution' + ': ' + str(1.0) + '\n'
        print(toPrint)

    def exit(self):
        """
        Des.
        """
        self.SU2.Output(1) # Temporary hack
        self.SU2.Postprocessing()
        print("***************************** Exit SU2 Solid solver *****************************")

class SU2SolidAdjoint(SU2SolidSolver, SolidAdjointSolver):
    """
    SU2 adjoint solver interface.
    """
    def initializeSolver(self, confFile, have_MPI, MPIComm):
        # --- Instantiate the single zone driver of SU2 --- #
        try:
            self.SU2 = pysu2.CDiscAdjSinglezoneDriver(confFile, 1, MPIComm)
        except TypeError as exception:
            print(('A TypeError occured in pysu2.CDiscAdjSinglezoneDriver : ',exception))
            if have_MPI == True:
                print('ERROR : You are trying to initialize MPI with a serial build of the wrapper. Please, remove the --parallel option that is incompatible with a serial build.')
            else:
                print('ERROR : You are trying to launch a computation without initializing MPI but the wrapper has been built in parallel. Please add the --parallel option in order to initialize MPI for the wrapper.')
    
    def initializeVariables(self):
        """
        Initialize variables required by the solver
        """
        SolidAdjointSolver.__init__(self)

    def applyNodalAdjointDisplacement(self, disp_adj_X, disp_adj_Y, disp_adj_Z, haloNodesDisplacements, time):
        PhysicalIndex = 0
        for iVertex in range(self.nNodes):
            GlobalIndex = self.SU2.GetVertexGlobalIndex(self.solidInterfaceID, iVertex)
            # in case of halo node, use haloNodesDisplacements with global fluid indexing
            if not GlobalIndex in list(self.haloNodeList.keys()):
                dispX = disp_adj_X[PhysicalIndex]
                dispY = disp_adj_Y[PhysicalIndex]
                dispZ = disp_adj_Z[PhysicalIndex]
                PhysicalIndex += 1
                self.SU2.SetSourceTerm_DispAdjoint(self.solidInterfaceID, iVertex, dispX, dispY, dispZ)

    def run(self, t1, t2):
        """
        Run one computation of SU2.
        """
        self.__steadyRun()

        self.__setCurrentState()

    def __setCurrentState(self):
        """
        Get the nodal (physical) displacements and velocities from SU2 structural solver.
        """

        PhysicalIndex = 0
        for iVertex in range(self.nNodes):
            # identify the halo nodes and ignore their nodal loads
            halo = self.SU2.ComputeVertexForces(self.solidInterfaceID, iVertex)
            self.SU2.ComputeVertexHeatFluxes(self.solidInterfaceID, iVertex)
            if halo == False:

                disp = self.SU2.GetFEA_Displacements(self.solidInterfaceID, iVertex)
                vel = self.SU2.GetFEA_Velocity(self.solidInterfaceID, iVertex)
                vel_n = self.SU2.GetFEA_Velocity_n(self.solidInterfaceID, iVertex)

                self.nodalDisp_X[PhysicalIndex] = disp[0]
                self.nodalDisp_Y[PhysicalIndex] = disp[1]
                self.nodalDisp_Z[PhysicalIndex] = disp[2]

                self.nodalVel_X[PhysicalIndex] = vel[0]
                self.nodalVel_Y[PhysicalIndex] = vel[1]
                self.nodalVel_Z[PhysicalIndex] = vel[2]

                self.nodalVel_XNm1[PhysicalIndex] = vel_n[0]
                self.nodalVel_YNm1[PhysicalIndex] = vel_n[1]
                self.nodalVel_ZNm1[PhysicalIndex] = vel_n[2]

                load = self.SU2.GetFlowLoad_Sensitivity(self.solidInterfaceID, iVertex)

                self.nodalAdjLoad_X[PhysicalIndex] = load[0]
                self.nodalAdjLoad_Y[PhysicalIndex] = load[1]
                self.nodalAdjLoad_Z[PhysicalIndex] = load[2]

                PhysicalIndex += 1

    def __steadyRun(self):
        """
        Run SU2 up to a converged steady state.
        """
        self.SU2.ResetConvergence()
        self.SU2.Preprocess(0)
        self.SU2.Run()
        StopIntegration = self.SU2.Monitor(0)
        self.SU2.Postprocess()
        self.SU2.Update()
        self.SU2.Output(0)
