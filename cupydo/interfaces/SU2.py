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

SU2.py
Python interface between the wrapper of SU2 and CUPyDO.
Authors D. THOMAS

'''

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

try: # Try to import the AD-capable SU2 module first
    import pysu2ad as pysu2
    adjoint = True
except ModuleNotFoundError as error:
    import pysu2
    adjoint = False
import numpy as np
from ..genericSolvers import FluidSolver, FluidAdjointSolver

# ----------------------------------------------------------------------
#  SU2 solver interface class
# ----------------------------------------------------------------------

class SU2(FluidSolver):
    """
    SU2 solver interface.
    """

    def __init__(self, p, have_MPI, MPIComm):
        """
        Initialize the SU2 solver and all the required interface variables.
        """

        self.initializeSolver(p['cfdFile'], have_MPI, MPIComm)
        allMovingMarkersTags = self.SU2.GetAllDeformMeshMarkersTag()                    # list containing the tags of all moving markers
        allCHTMarkersTags = self.SU2.GetAllCHTMarkersTag()
        allMarkersID = self.SU2.GetAllBoundaryMarkers()                             # dic : allMarkersID['marker_tag'] = marker_ID
        self.fluidInterfaceID = None                                                # identification of the f/s boundary, currently limited to one boundary, by default the first tag in allMovingMarkersTags
        if not allMovingMarkersTags and not allCHTMarkersTags:
            #raise Exception('No interface for FSI was defined.')
            self.fluidInterfaceID = None
        elif allMovingMarkersTags and not allCHTMarkersTags:
            if allMovingMarkersTags[0] in list(allMarkersID.keys()):
                self.fluidInterfaceID = allMarkersID[allMovingMarkersTags[0]]
        elif not allMovingMarkersTags and allCHTMarkersTags:
            if allCHTMarkersTags[0] in list(allMarkersID.keys()):
                self.fluidInterfaceID = allMarkersID[allCHTMarkersTags[0]]
        elif allMovingMarkersTags and allCHTMarkersTags:
            if allMovingMarkersTags[0] == allCHTMarkersTags[0]:
                if allMovingMarkersTags[0] in list(allMarkersID.keys()):
                    self.fluidInterfaceID = allMarkersID[allMovingMarkersTags[0]]
            else:
                raise Exception("Moving and CHT markers have to be the same!\n")

        self.regime = p['regime']                                   # computation type : steady or unsteady

        # --- Calculate the number of nodes (on each partition) --- #
        self.nNodes = 0
        self.nHaloNode = 0
        self.nPhysicalNodes = 0
        if self.fluidInterfaceID != None:
            self.nNodes = self.SU2.GetNumberVertices(self.fluidInterfaceID)           # numbers of nodes at the f/s interface (halo+physical)
            self.nHaloNode = self.SU2.GetNumberHaloVertices(self.fluidInterfaceID)    # numbers of nodes at the f/s interface (halo)
        self.nPhysicalNodes = self.nNodes - self.nHaloNode                        # numbers of nodes at the f/s interface (physical)
        # Initialize list to store point indices
        self.pointIndexList = np.zeros(self.nPhysicalNodes, dtype=int)

        # --- initial position of the f/s interface(s) --- #
        self.nodalInitialPos_X = np.zeros((self.nPhysicalNodes))             # initial position of the f/s interface
        self.nodalInitialPos_Y = np.zeros((self.nPhysicalNodes))
        self.nodalInitialPos_Z = np.zeros((self.nPhysicalNodes))
        self.haloNodesPositionsInit = {}

        self.initializeVariables(p)

        # --- Initialize the interface position and the nodal loads --- #
        PhysicalIndex = 0
        for iVertex in range(self.nNodes):
            posX, posY, posZ = self.SU2.GetInitialMeshCoord(self.fluidInterfaceID, iVertex)
            # posX = self.SU2.GetVertexCoordX(self.fluidInterfaceID, iVertex)
            # posY = self.SU2.GetVertexCoordY(self.fluidInterfaceID, iVertex)
            # posZ = self.SU2.GetVertexCoordZ(self.fluidInterfaceID, iVertex)
            if self.SU2.IsAHaloNode(self.fluidInterfaceID, iVertex):
                GlobalIndex = self.SU2.GetVertexGlobalIndex(self.fluidInterfaceID, iVertex)
                self.haloNodeList[GlobalIndex] = iVertex
                self.haloNodesPositionsInit[GlobalIndex] = (posX, posY, posZ)
            else:
                GlobalIndex = self.SU2.GetVertexGlobalIndex(self.fluidInterfaceID, iVertex)
                self.pointIndexList[PhysicalIndex] = GlobalIndex
                # self.SU2.ComputeVertexForces(self.fluidInterfaceID, iVertex)
                Fx, Fy, Fz = self.SU2.GetFlowLoad(self.fluidInterfaceID, iVertex)
                Temp = self.SU2.GetVertexTemperature(self.fluidInterfaceID, iVertex)
                self.nodalInitialPos_X[PhysicalIndex] = posX
                self.nodalInitialPos_Y[PhysicalIndex] = posY
                self.nodalInitialPos_Z[PhysicalIndex] = posZ
                self.nodalLoad_X[PhysicalIndex] = Fx
                self.nodalLoad_Y[PhysicalIndex] = Fy
                self.nodalLoad_Z[PhysicalIndex] = Fz
                self.nodalTemperature[PhysicalIndex] = Temp
                PhysicalIndex += 1

        self.initRealTimeData()

        # print("\n -------------------------- FLUID NODES ------------------------------ \n")
        # print("There is a total of", self.nNodes, "nodes\n")
        # for iIndex in range(self.nPhysicalNodes):
        #     print(self.pointIndexList[iIndex], self.nodalInitialPos_X[iIndex], self.nodalInitialPos_Y[iIndex], self.nodalInitialPos_Z[iIndex])

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

    def initializeVariables(self, p):
        """
        Initialize variables required by the solver
        """
        FluidSolver.__init__(self, p)

    def run(self, t1, t2):
        """
        Run one computation of SU2.
        """

        if self.regime == 'unsteady':

            dt = t2-t1
            if not np.allclose(self.SU2.GetUnsteady_TimeStep(),dt):
                raise Exception('SU2 and FSI time step do not match')
            nt = int(t2/dt)
            self.__unsteadyRun(nt)
        else:
            self.__steadyRun()

        self.__setCurrentState()
        return True

    def __unsteadyRun(self, nt):
        """
        Run SU2 on one time step.
        """

        self.SU2.ResetConvergence()
        self.SU2.Run()
        self.SU2.Postprocess()

    def __steadyRun(self):
        """
        Run SU2 up to a converged steady state.
        """

        self.SU2.ResetConvergence()
        self.SU2.Preprocess(0)
        self.SU2.Run()
        StopIntegration = self.SU2.Monitor(0)
        self.SU2.Postprocess()
        self.SU2.Output(0)

    def __setCurrentState(self):
        """
        Get the nodal (physical) loads from SU2 solver.
        """

        PhysicalIndex = 0
        for iVertex in range(self.nNodes):
            # identify the halo nodes and ignore their nodal loads
            # halo = self.SU2.ComputeVertexForces(self.fluidInterfaceID, iVertex)
            halo = self.SU2.IsAHaloNode(self.fluidInterfaceID, iVertex)
            # self.SU2.ComputeVertexHeatFluxes(self.fluidInterfaceID, iVertex)
            if halo == False:
                Fx, Fy, Fz = self.SU2.GetFlowLoad(self.fluidInterfaceID, iVertex)
                Temp = self.SU2.GetVertexTemperature(self.fluidInterfaceID, iVertex)
                WallHF = self.SU2.GetVertexNormalHeatFlux(self.fluidInterfaceID, iVertex)
                Qx, Qy, Qz = self.SU2.GetVertexHeatFluxes(self.fluidInterfaceID, iVertex)
                self.nodalLoad_X[PhysicalIndex] = Fx
                self.nodalLoad_Y[PhysicalIndex] = Fy
                self.nodalLoad_Z[PhysicalIndex] = Fz
                self.nodalTemperature[PhysicalIndex] = Temp
                self.nodalNormalHeatFlux[PhysicalIndex] = WallHF
                self.nodalHeatFlux_X[PhysicalIndex] = Qx
                self.nodalHeatFlux_Y[PhysicalIndex] = Qy
                self.nodalHeatFlux_Z[PhysicalIndex] = Qz
                PhysicalIndex += 1

    def getNodalIndex(self, iVertex):
        """
        Return the global index (fluid solver index) of a node.
        """

        return self.SU2.GetVertexGlobalIndex(self.fluidInterfaceID, iVertex)

    def getNodalInitialPositions(self):
        """
        Return the initial position of the f/s interface.
        """

        return (self.nodalInitialPos_X, self.nodalInitialPos_Y, self.nodalInitialPos_Z)

    def applyNodalDisplacements(self, disp_X, disp_Y, disp_Z, dt, haloNodesDisplacements):
        """
        Set the displacement of the f/s boundary before mesh morphing.
        """

        PhysicalIndex = 0
        for iVertex in range(self.nNodes):
            GlobalIndex = self.SU2.GetVertexGlobalIndex(self.fluidInterfaceID, iVertex)
            # in case of halo node, use haloNodesDisplacements with global fluid indexing
            if GlobalIndex in list(self.haloNodeList.keys()):
                dispX, dispY, dispZ = haloNodesDisplacements[GlobalIndex]
                posX0, posY0, posZ0 = self.haloNodesPositionsInit[GlobalIndex]
                newPosX = posX0 + dispX
                newPosY = posY0 + dispY
                newPosZ = posZ0 + dispZ
            else:
                dispX = disp_X[PhysicalIndex]
                dispY = disp_Y[PhysicalIndex]
                dispZ = disp_Z[PhysicalIndex]
                newPosX = dispX + self.nodalInitialPos_X[PhysicalIndex]
                newPosY = dispY + self.nodalInitialPos_Y[PhysicalIndex]
                newPosZ = dispZ + self.nodalInitialPos_Z[PhysicalIndex]
                PhysicalIndex += 1
            self.SU2.SetMeshDisplacement(self.fluidInterfaceID, iVertex, dispX, dispY, dispZ)

    def applyNodalHeatFluxes(self, HF_X, HF_Y, HF_Z, dt, haloNodesHeatFlux):
        """
        Set the heat fluxes on the f/s boundary and update the multi-grid structure (if any).
        """

        PhysicalIndex = 0
        for iVertex in range(self.nNodes):
            WallHF = 0.0
            if not self.SU2.IsAHaloNode(self.fluidInterfaceID, iVertex):
                N = self.SU2.GetVertexUnitNormal(self.fluidInterfaceID, iVertex)
                #In SU2, the surface normal is pointing inwards the fluid domain (meaning outwards the solid domain).
                WallHF = HF_X[PhysicalIndex]*N[0] + HF_Y[PhysicalIndex]*N[1] + HF_Z[PhysicalIndex]*N[2]
                self.SU2.SetVertexNormalHeatFlux(self.fluidInterfaceID, iVertex, WallHF)
                PhysicalIndex += 1

    def applyNodalTemperatures(self, Temperature, dt, haloNodesTemperature):


        PhysicalIndex = 0
        for iVertex in range(self.nNodes):
            if not self.SU2.IsAHaloNode(self.fluidInterfaceID, iVertex):
                self.SU2.SetVertexTemperature(self.fluidInterfaceID, iVertex, Temperature[PhysicalIndex])
                PhysicalIndex += 1

    def update(self, dt):
        """
        Update the solution after each time step (converged) for an unsteady computation.
        """
        
        self.SU2.Update()

    def save(self, nt):
        """
        Write the solution on disk.
        """

        stopComp = self.SU2.Monitor(nt)
        self.SU2.Output(nt)

        return stopComp

    def initRealTimeData(self):

        solFile = open('AerodynamicCoeff.ascii', "w")
        solFile.write("{0:>12s}   {1:>12s}   {2:>12s}   {3:>12s}\n".format("Time", "FSI_Iter", "C_Lift", "C_Drag"))
        solFile.close()
    
    def saveRealTimeData(self, time, nFSIIter):

        CLift = self.SU2.Get_LiftCoeff()
        CDrag = self.SU2.Get_DragCoeff()

        solFile = open('AerodynamicCoeff.ascii', "a")
        solFile.write("{0:12.6f}   {1:12d}   {2:12.6f}   {3:12.6f}\n".format(time, nFSIIter, CLift, CDrag))
        solFile.close()

    def printRealTimeData(self, time, nFSIIter):


        CLift = self.SU2.Get_LiftCoeff()
        CDrag = self.SU2.Get_DragCoeff()

        print('RES-FSI-Lift coeff : '+ str(CLift))
        print('RES-FSI-Drag coeff : '+ str(CDrag))

    def meshUpdate(self, nt):
        """
        Perform the mesh morphing.
        """

        if self.regime == 'unsteady' and nt>0:
            self.SU2.DynamicMeshUpdate(nt)
        elif self.regime == 'steady':
            self.SU2.StaticMeshUpdate()

    def boundaryConditionsUpdate(self):

        self.SU2.BoundaryConditionsUpdate()

    def setInitialMeshDeformation(self):
        """
        Perform the mesh morphing for FSI initial conditions.
        """

        self.SU2.StaticMeshUpdate()

    def preprocessTimeIter(self, timeIter):
        """
        Preprocessing routine before each time step.
        """

        self.SU2.Preprocess(timeIter)
 
    def exit(self):
        """
        Exit and clean the SU2 solver.
        """

        self.SU2.Postprocessing()

    def fakeSolidSolver(self, dt):

        return

class SU2Adjoint(SU2, FluidAdjointSolver):
    """
    SU2 adjoint solver interface.
    """
    def initializeSolver(self, confFile, have_MPI, MPIComm):
        # --- Instantiate the single zone driver of SU2 --- #
        if not adjoint:
            print('ERROR: You are trying to launch an adjoint calculation with an AD-incapable build of the SU2 wrapper. Please, add the -Denable-pywrapper=true to meson options.')
            return
        try:
            self.SU2 = pysu2.CDiscAdjSinglezoneDriver(confFile, 1, MPIComm)
        except TypeError as exception:
            print(('A TypeError occured in pysu2.CDiscAdjSinglezoneDriver : ',exception))
            if have_MPI == True:
                print('ERROR : You are trying to initialize MPI with a serial build of the wrapper. Please, remove the --parallel option that is incompatible with a serial build.')
            else:
                print('ERROR : You are trying to launch a computation without initializing MPI but the wrapper has been built in parallel. Please add the --parallel option in order to initialize MPI for the wrapper.')
    
    def initializeVariables(self, p):
        """
        Initialize variables required by the solver
        """
        FluidAdjointSolver.__init__(self, p)

    def __setCurrentState(self):
        SU2._SU2__setCurrentState(self) # Using SU2 original __setCurrentState
        PhysicalIndex = 0
        for iVertex in range(self.nNodes):
            # identify the halo nodes and ignore their nodal loads
            # halo = self.SU2.ComputeVertexForces(self.fluidInterfaceID, iVertex)
            halo = self.SU2.IsAHaloNode(self.fluidInterfaceID, iVertex)
            if halo == False:
                self.nodalAdjDisp_X[PhysicalIndex], self.nodalAdjDisp_Y[PhysicalIndex], self.nodalAdjDisp_Z[PhysicalIndex]  = self.SU2.GetMeshDisp_Sensitivity(self.fluidInterfaceID, iVertex)
                PhysicalIndex += 1

    def run(self, t1, t2):
        """
        Run one computation of SU2.
        """
        self.__steadyRun()

        self.__setCurrentState()
        return True

    def __steadyRun(self):
        """
        Run SU2 up to a converged steady state.
        """

        self.SU2.ResetConvergence()
        self.SU2.Preprocess(0)
        self.SU2.Run()
        self.SU2.Postprocess()
        self.SU2.Update()
        StopIntegration = self.SU2.Monitor(0)
        self.SU2.Output(0)
    
    def applyNodalAdjointForce(self, load_X, load_Y, load_Z, dt, haloNodesLoads):
        PhysicalIndex = 0
        for iVertex in range(self.nNodes):
            GlobalIndex = self.SU2.GetVertexGlobalIndex(self.fluidInterfaceID, iVertex)
            # in case of halo node, use haloNodesLoads with global fluid indexing
            if GlobalIndex in list(self.haloNodeList.keys()):
                loadX, loadY, loadZ = haloNodesLoads[GlobalIndex]
            else:
                loadX = load_X[PhysicalIndex]
                loadY = load_Y[PhysicalIndex]
                loadZ = load_Z[PhysicalIndex]
                PhysicalIndex += 1
            self.SU2.SetFlowLoad_Adjoint(self.fluidInterfaceID, iVertex, loadX, loadY, loadZ)
