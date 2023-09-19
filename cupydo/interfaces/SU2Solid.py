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

SU2Solid.py
Python interface between the wrapper of SU2 for solid mechanics and CUPyDO.
Authors R. Sanchez - TU Kaiserslautern

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
import math
import numpy as np
try: # Try to import the optional vtk module. If unavailable, send a warning
    import vtk
except ModuleNotFoundError as error:
    vtk = None
    UserWarning('vtk module not found! SU2Solid solution output will not work\n')

# Those are mandatory
import numpy as np
from ..genericSolvers import SolidSolver, SolidAdjointSolver


# ----------------------------------------------------------------------
#  SU2SolidSolver class
# ----------------------------------------------------------------------

class SU2SolidSolver(SolidSolver):
    """
    SU2 solver interface.
    """

    def __init__(self, p, have_MPI, MPIComm=None):
        """
        Initialize the SU2 solver and all the required interface variables.
        """

        self.initializeSolver(p['csdFile'], have_MPI, MPIComm)
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

        self.computationType = p['compType']  # computation type : steady (default) or unsteady
        self.nodalLoadsType = p['nodalLoadsType']  # nodal loads type to extract : force (in N, default) or pressure (in Pa)

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

        self.extractors = p['extractors'] # List of points to extract 
        self.filename = p['surfaceFilename'] # Filename (no extension) of surface output file
        self.extension = p['surfaceExtension'] # Extension of surface output file

        self.__setCurrentState()
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

            dt = t2-t1
            if not np.allclose(self.SU2.GetUnsteady_TimeStep(),dt):
                raise Exception('SU2 and FSI time step do not match')
            self.__unsteadyRun(t1, t2)
        else:
            self.__steadyRun()

        self.__setCurrentState()
        return True

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

                PhysicalIndex += 1

    def getNodalInitialPositions(self):


        return (self.nodalInitialPos_X, self.nodalInitialPos_Y, self.nodalInitialPos_Z)

    def getNodalIndex(self, iVertex):
        """
        Returns the index (identifier) of the iVertex^th interface node.
        """

        return self.SU2.GetVertexGlobalIndex(self.solidInterfaceID, int(iVertex))

    def getNodalDisplacements(self):


        return (self.nodalDisp_X, self.nodalDisp_Y, self.nodalDisp_Z)

    def applyNodalForce(self, load_X, load_Y, load_Z, dt, haloNodesLoads = {}):


        # --- Initialize the interface position and the nodal loads --- #
        PhysicalIndex = 0
        for iVertex in range(self.nNodes):
            GlobalIndex = self.SU2.GetVertexGlobalIndex(self.solidInterfaceID, iVertex)
            # in case of halo node, use haloNodesDisplacements with global fluid indexing
            if GlobalIndex in list(self.haloNodeList.keys()):
                LoadX, LoadY, LoadZ = haloNodesLoads[GlobalIndex]
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

    def save(self):


        stopComp = self.SU2.Monitor(1)
        self.SU2.Output(1)

        return stopComp

    def initRealTimeData(self):


        if not vtk is None:
            if self.extension == 'vtu':
                self.dataReader = vtk.vtkXMLUnstructuredGridReader()
                self.dataReader.SetFileName(self.filename + '.' + self.extension)
            elif self.extension == 'vtk':
                self.dataReader = vtk.vtkUnstructuredGridReader()
                self.dataReader.ReadAllScalarsOn()
                self.dataReader.ReadAllVectorsOn()
                self.dataReader.ReadAllTensorsOn()
                self.dataReader.ReadAllFieldsOn()
                self.dataReader.SetFileName(self.filename + '.' + self.extension)
            elif self.extension == 'dat':
                self.dataReader = vtk.vtkTecplotReader()
                self.dataReader.SetFileName(self.filename + '.' + self.extension)
            else:
                UserWarning('Reader for extension {} not available. SU2Solid solution output will not work\n'.format(self.extension))
                self.dataReader = None
        else:
            self.dataReader = None
        
        solFile = open('SolidSolution.ascii', "w")
        solFile.write('{:>12s}   {:>12s}'.format('Time', 'Iteration'))
        for gidx in self.extractors:
            solFile.write('   {:>12s}   {:>12s}   {:>12s}'.format('x_'+str(gidx), 'y_'+str(gidx), 'z_'+str(gidx)))
        solFile.write('\n')
        solFile.close()

    def saveRealTimeData(self, time, nFSIIter):


        solFile = open('SolidSolution.ascii', "a")
        solFile.write("{:>12.6f}   {:>12d}".format(time, nFSIIter))
        if not self.dataReader is None:
            self.dataReader.Modified() # Update the dataReader with the new VTK file
            self.dataReader.Update()
            nodalDisp = np.array(self.dataReader.GetOutput().GetPointData().GetArray('Displacement'))
        else:
            nodalDisp = np.zeros((self.nNodes, 3)) # If no vtk module use zeros
        for gidx in self.extractors:
            solFile.write('   {:>12.10f}   {:>12.10f}   {:>12.10f}'.format(nodalDisp[gidx, 0], nodalDisp[gidx, 1], nodalDisp[gidx, 2]))
        solFile.write('\n')
        solFile.close()

    def printRealTimeData(self, time, nFSIIter):


        toPrint = 'RES-FSI-' + 'SU2SolidSolution' + ': ' + str(1.0) + '\n'
        print(toPrint)

    def exit(self):

        self.SU2.Output(1) # Temporary hack
        self.SU2.Postprocessing()
        print("***************************** Exit SU2 Solid solver *****************************")

class SU2SolidAdjoint(SU2SolidSolver, SolidAdjointSolver):
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
    
    def initializeVariables(self):
        """
        Initialize variables required by the solver
        """
        SolidAdjointSolver.__init__(self)

    def applyNodalAdjointDisplacement(self, disp_adj_X, disp_adj_Y, disp_adj_Z, haloNodesDisplacements, dt):
        PhysicalIndex = 0
        for iVertex in range(self.nNodes):
            GlobalIndex = self.SU2.GetVertexGlobalIndex(self.solidInterfaceID, iVertex)
            # in case of halo node, use haloNodesDisplacements with global fluid indexing
            if GlobalIndex in list(self.haloNodeList.keys()):
                dispX, dispY, dispZ = haloNodesDisplacements[GlobalIndex]
            else:
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
        return True

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
