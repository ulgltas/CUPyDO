#!/usr/bin/env python
# -*- coding: latin-1; -*-

# \file SU2Interface.py
#  \brief Python interface between the wrapper of SU2 and the FSI coupler.
#  \author D. THOMAS, University of Liege, Belgium. Department of Aerospace and Mechanical Engineering.
#  \version BETA

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import pysu2
import math
import numpy as np
from FSICoupler import FluidSolver

# ----------------------------------------------------------------------
#  SU2 solver class
# ----------------------------------------------------------------------

class SU2Solver(FluidSolver):
    """
    SU2 solver interface.
    """

    def __init__(self, confFile, bndno, nDim, computationType, nodalLoadsType, have_MPI, MPIComm=None):
        """
        Initialize the SU2 solver and all the required interface variables.
        """

        # --- Instantiate the single zone driver of SU2 --- #
        try:
            self.SU2 = pysu2.CFluidDriver(confFile, 1, nDim, MPIComm)
        except TypeError as exception:
            print('A TypeError occured in pysu2.CSingleZoneDriver : ',exception)
            if have_MPI == True:
                print('ERROR : You are trying to initialize MPI with a serial build of the wrapper. Please, remove the --parallel option that is incompatible with a serial build.')
            else:
                print('ERROR : You are trying to launch a computation without initializing MPI but the wrapper has been built in parallel. Please add the --parallel option in order to initialize MPI for the wrapper.')

        allMovingMarkersTags = self.SU2.GetAllMovingMarkersTag()                    # list containing the tags of all moving markers
        allMarkersID = self.SU2.GetAllBoundaryMarkers()                             # dic : allMarkersID['marker_tag'] = marker_ID
        self.fluidInterfaceID = None
        if allMovingMarkersTags[0] in allMarkersID.keys():
            self.fluidInterfaceID = allMarkersID[allMovingMarkersTags[0]]         # identification of the f/s boundary, currently limited to one boundary, by default the first tag in allMovingMarkersTags
        self.computationType = computationType                                    # computation type : steady (default) or unsteady
        self.nodalLoadsType = nodalLoadsType                                      # nodal loads type to extract : force (in N, default) or pressure (in Pa)

        # --- Calculate the number of nodes (on each partition) --- #
        self.nNodes = 0
        self.nHaloNode = 0
        self.nPhysicalNodes = 0
        if self.fluidInterfaceID != None:
            self.nNodes = self.SU2.GetNumberVertices(self.fluidInterfaceID)           # numbers of nodes at the f/s interface (halo+physical)
            self.nHaloNode = self.SU2.GetNumberHaloVertices(self.fluidInterfaceID)    # numbers of nodes at the f/s interface (halo)
        self.nPhysicalNodes = self.nNodes - self.nHaloNode                        # numbers of nodes at the f/s interface (physical)
    
        self.nodalInitialPos_X = np.zeros((self.nPhysicalNodes))             # initial position of the f/s interface
        self.nodalInitialPos_Y = np.zeros((self.nPhysicalNodes))
        self.nodalInitialPos_Z = np.zeros((self.nPhysicalNodes))
        self.haloNodesPositionsInit = {}

        FluidSolver.__init__(self)

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
                self.SU2.ComputeVertexForces(self.fluidInterfaceID, iVertex)
                Fx = self.SU2.GetVertexForceX(self.fluidInterfaceID, iVertex)
                Fy = self.SU2.GetVertexForceY(self.fluidInterfaceID, iVertex)
                Fz = self.SU2.GetVertexForceZ(self.fluidInterfaceID, iVertex)
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
        NbIter = self.SU2.GetnExtIter()
        Iter = 0
        while Iter < NbIter:
            self.SU2.PreprocessExtIter(Iter)
            self.SU2.Run()
            StopIntegration = self.SU2.Monitor(Iter)
            self.SU2.Output(Iter)
            if StopIntegration:
                break;
            Iter += 1

    def __setCurrentState(self):
        """
        Get the nodal (physical) loads from SU2 solver.
        """

        PhysicalIndex = 0
        for iVertex in range(self.nNodes):
            # identify the halo nodes and ignore their nodal loads
            halo = self.SU2.ComputeVertexForces(self.fluidInterfaceID, iVertex)
            self.SU2.ComputeVertexHeatFluxes(self.fluidInterfaceID, iVertex)
            if halo == False:
                if self.nodalLoadsType == 'pressure':
                    Fx = self.SU2.GetVertexForceDensityX(self.fluidInterfaceID, iVertex)
                    Fy = self.SU2.GetVertexForceDensityY(self.fluidInterfaceID, iVertex)
                    Fz = self.SU2.GetVertexForceDensityZ(self.fluidInterfaceID, iVertex)
                else:
                    Fx = self.SU2.GetVertexForceX(self.fluidInterfaceID, iVertex)
                    Fy = self.SU2.GetVertexForceY(self.fluidInterfaceID, iVertex)
                    Fz = self.SU2.GetVertexForceZ(self.fluidInterfaceID, iVertex)
                Temp = self.SU2.GetVertexTemperature(self.fluidInterfaceID, iVertex)
                WallHF = self.SU2.GetVertexNormalHeatFlux(self.fluidInterfaceID, iVertex)
                Qx = self.SU2.GetVertexHeatFluxX(self.fluidInterfaceID, iVertex)
                Qy = self.SU2.GetVertexHeatFluxY(self.fluidInterfaceID, iVertex)
                Qz = self.SU2.GetVertexHeatFluxZ(self.fluidInterfaceID, iVertex)
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

    def applyNodalDisplacements(self, disp_X, disp_Y, disp_Z, dispnM1_X, dispnM1_Y, dispnM1_Z, haloNodesDisplacements, time):
        """
        Set the displacement of the f/s boundary before mesh morphing.
        """

        PhysicalIndex = 0
        for iVertex in range(self.nNodes):
            GlobalIndex = self.SU2.GetVertexGlobalIndex(self.fluidInterfaceID, iVertex)
            # in case of halo node, use haloNodesDisplacements with global fluid indexing
            if GlobalIndex in self.haloNodeList.keys():
                dispX, dispY, dispZ = haloNodesDisplacements[GlobalIndex]
                posX0, posY0, posZ0 = self.haloNodesPositionsInit[GlobalIndex]
                newPosX = posX0 + dispX
                newPosY = posY0 + dispY
                newPosZ = posZ0 + dispZ
            else:
                newPosX = disp_X[PhysicalIndex] + self.nodalInitialPos_X[PhysicalIndex]
                newPosY = disp_Y[PhysicalIndex] + self.nodalInitialPos_Y[PhysicalIndex]
                newPosZ = disp_Z[PhysicalIndex] + self.nodalInitialPos_Z[PhysicalIndex]
                PhysicalIndex += 1
            self.SU2.SetVertexCoordX(self.fluidInterfaceID, iVertex, newPosX)
            self.SU2.SetVertexCoordY(self.fluidInterfaceID, iVertex, newPosY)
            self.SU2.SetVertexCoordZ(self.fluidInterfaceID, iVertex, newPosZ)
            self.SU2.SetVertexVarCoord(self.fluidInterfaceID, iVertex)

    def applyNodalHeatFluxes(self, HF_X, HF_Y, HF_Z, time):
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

        if self.fluidInterfaceID != None:
            self.SU2.MGUpdateBoundaryConditions_HeatFlux(self.fluidInterfaceID)

    def applyNodalTemperatures(self, Temperature, time):
        """
        Des.
        """

        PhysicalIndex = 0
        for iVertex in range(self.nNodes):
            if not self.SU2.IsAHaloNode(self.fluidInterfaceID, iVertex):
                self.SU2.SetVertexTemperature(self.fluidInterfaceID, iVertex, Temperature[PhysicalIndex])
                PhysicalIndex += 1

        if self.fluidInterfaceID != None:
            self.SU2.MGUpdateBoundaryConditions_Temperature(self.fluidInterfaceID)

    def setInitialInterfaceHeatFlux(self):
        """
        Set an initial (first guess) and uniform heat flux on the f/s boundary.
        """

        PhysicalIndex = 0
        for iVertex in range(self.nNodes):
            if not self.SU2.IsAHaloNode(self.fluidInterfaceID, iVertex):
                self.SU2.SetVertexNormalHeatFlux(self.fluidInterfaceID, iVertex, self.QWallInit)
                PhysicalIndex += 1

        if self.fluidInterfaceID != None:
            self.SU2.MGUpdateBoundaryConditions_HeatFlux(self.fluidInterfaceID)

    def setInitialInterfaceTemperature(self):
        """
        Des
        """

        PhysicalIndex = 0
        for iVertex in range(self.nNodes):
            if not self.SU2.IsAHaloNode(self.fluidInterfaceID, iVertex):
                self.SU2.SetVertexTemperature(self.fluidInterfaceID, iVertex, self.TWallInit)
                PhysicalIndex += 1

        if self.fluidInterfaceID != None:
            self.SU2.MGUpdateBoundaryConditions_Temperature(self.fluidInterfaceID)

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
        """
        Description.
        """

        solFile = open('AerodynamicCoeff.ascii', "w")
        solFile.write("Time\tnIter\tCl\tCd\n")
        solFile.close()
    
    def saveRealTimeData(self, time, nFSIIter):
        """
        Description.
        """

        CLift = self.SU2.Get_LiftCoeff()
        CDrag = self.SU2.Get_DragCoeff()

        solFile = open('AerodynamicCoeff.ascii', "a")
        solFile.write(str(time) + '\t' + str(nFSIIter) + '\t' + str(CLift)  + '\t' + str(CDrag) + '\n')

    def meshUpdate(self, nt):
        """
        Perform the mesh morphing.
        """

        if self.computationType == 'unsteady':
            self.SU2.DynamicMeshUpdate(nt)
        else:
            self.SU2.StaticMeshUpdate()

    def setInitialMeshDeformation(self):
        """
        Perform the mesh morphing for FSI initial conditions.
        """

        self.SU2.SetInitialMesh()

    def preprocessTimeIter(self, timeIter):
        """
        Preprocessing routine before each time step.
        """

        self.SU2.PreprocessExtIter(timeIter)

    def remeshing(self):
        """
        Desctiption.
        """

        return
 
    def exit(self):
        """
        Exit and clean the SU2 solver.
        """

        self.SU2.Postprocessing()

    def fakeSolidSolver(self, time):
        """
        Des.
        """

        return
