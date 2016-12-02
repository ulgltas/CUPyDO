#!/usr/bin/env python
# -*- coding: latin-1; -*-

# \file SU2Interface.py
#  \brief Description.
#  \author D. THOMAS, University of Liege, Belgium. Department of Aerospace and Mechanical Engineering
#  \version BETA

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import WrapSU2
import math
import numpy as np
from fsi import FluidSolver

# ----------------------------------------------------------------------
#  SU2 solver class
# ----------------------------------------------------------------------

class SU2Solver(FluidSolver):
  def __init__(self, confFile, bndno, nDim, computationType, have_MPI, MPIComm):
    """
    Description.
    """

    try:
      self.SU2 = WrapSU2.CSingleZoneDriver(confFile, 1, nDim, MPIComm)
    except TypeError as exception:
      print('A TypeError occured in WrapSU2.CSingleZoneDriver : ',exception)
      if have_MPI == True:
        print('ERROR : You are trying to initialize MPI with a serial build of the wrapper. Please, remove the --parallel option that is incompatible with a serial build.')
      else:
        print('ERROR : You are trying to launch a computation without initializing MPI but the wrapper has been built in parallel. Please add the --parallel option in order to initialize MPI for the wrapper.')

    self.fluidInterfaceID = self.SU2.GetMovingMarker()
    self.computationType = computationType

    # Calculate the number of physical (= not halo) nodes on each partition
    self.nNodes = self.SU2.GetNumberVertices(self.fluidInterfaceID)
    self.nHaloNode = self.SU2.GetNumberHaloVertices(self.fluidInterfaceID)
    self.nPhysicalNodes = self.nNodes - self.nHaloNode
    
    self.interface_array_X_init = np.zeros((self.nPhysicalNodes))
    self.interface_array_Y_init = np.zeros((self.nPhysicalNodes))
    self.interface_array_Z_init = np.zeros((self.nPhysicalNodes))
    self.haloNodesPositionsInit = {}

    FluidSolver.__init__(self)

    localIndex = 0
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
        self.interface_array_X_init[localIndex] = posX
        self.interface_array_Y_init[localIndex] = posY
        self.interface_array_Z_init[localIndex] = posZ
        self.nodalLoad_X[localIndex] = Fx
        self.nodalLoad_Y[localIndex] = Fy
        self.nodalLoad_Z[localIndex] = Fz
        localIndex += 1

  def run(self, t1, t2):

    if self.computationType == 'unsteady':
      self.__unsteadyRun(t1, t2)
    else:
      self.__steadyRun()

  def __unsteadyRun(self, t1, t2):
    """
    Description.
    """

    self.SU2.ResetConvergence()
    self.SU2.Run()

    localIndex = 0
    for iVertex in range(self.nNodes):
      halo = self.SU2.ComputeVertexForces(self.fluidInterfaceID, iVertex)
      if halo == False:
        Fx = self.SU2.GetVertexForceX(self.fluidInterfaceID, iVertex)
        Fy = self.SU2.GetVertexForceY(self.fluidInterfaceID, iVertex)
        Fz = self.SU2.GetVertexForceZ(self.fluidInterfaceID, iVertex)
        self.nodalLoad_X[localIndex] = Fx
        self.nodalLoad_Y[localIndex] = Fy
        self.nodalLoad_Z[localIndex] = Fz
        localIndex += 1

  def __steadyRun(self):
    """
    Descr.
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

    localIndex = 0
    for iVertex in range(self.nNodes):
      halo = self.SU2.ComputeVertexForces(self.fluidInterfaceID, iVertex)
      if halo == False:
        Fx = self.SU2.GetVertexForceX(self.fluidInterfaceID, iVertex)
        Fy = self.SU2.GetVertexForceY(self.fluidInterfaceID, iVertex)
        Fz = self.SU2.GetVertexForceZ(self.fluidInterfaceID, iVertex)
        self.nodalLoad_X[localIndex] = Fx
        self.nodalLoad_Y[localIndex] = Fy
        self.nodalLoad_Z[localIndex] = Fz
        localIndex += 1

  def getNodalIndex(self, iVertex):
    """
    Description.
    """

    return self.SU2.GetVertexGlobalIndex(self.fluidInterfaceID, iVertex)

  def getNodalInitialPositions(self):
    """
    Description.
    """

    return (self.interface_array_X_init, self.interface_array_Y_init, self.interface_array_Z_init)

  def applyNodalDisplacements(self, disp_X, disp_Y, disp_Z, dispnM1_X, dispnM1_Y, dispnM1_Z, haloNodesDisplacements, time):
    """
    Description
    """

    PhysicalIndex = 0
    for iVertex in range(self.nNodes):
      GlobalIndex = self.SU2.GetVertexGlobalIndex(self.fluidInterfaceID, iVertex)
      if GlobalIndex in self.haloNodeList.keys():
        dispX, dispY, dispZ = haloNodesDisplacements[GlobalIndex]
        posX0, posY0, posZ0 = self.haloNodesPositionsInit[GlobalIndex]
        newPosX = posX0 + dispX
        newPosY = posY0 + dispY
        newPosZ = posZ0 + dispZ
      else:
        newPosX = disp_X[PhysicalIndex] + self.interface_array_X_init[PhysicalIndex]
        newPosY = disp_Y[PhysicalIndex] + self.interface_array_Y_init[PhysicalIndex]
        newPosZ = disp_Z[PhysicalIndex] + self.interface_array_Z_init[PhysicalIndex]
        PhysicalIndex += 1
      self.SU2.SetVertexCoordX(self.fluidInterfaceID, iVertex, newPosX)
      self.SU2.SetVertexCoordY(self.fluidInterfaceID, iVertex, newPosY)
      self.SU2.SetVertexCoordZ(self.fluidInterfaceID, iVertex, newPosZ)
      nodalVarCoordNorm = self.SU2.SetVertexVarCoord(self.fluidInterfaceID, iVertex)


  def update(self, dt):
    """
    Description.
    """

    self.SU2.Update()

  def save(self, nt):
    """
    Description.
    """

    stopComp = self.SU2.Monitor(nt)
    self.SU2.Output(nt)

    return stopComp

    def initRealTimeData(self):
        return
    
    def saveRealTimeData(self, time, nFSIIter):
        return

  def meshUpdate(self, nt):
    """
    Description.
    """

    if self.computationType == 'unsteady':
      self.SU2.DynamicMeshUpdate(nt)
    else:
      self.SU2.StaticMeshUpdate()

  def setInitialMeshDeformation(self):
    """
    Des.
    """

    self.SU2.SetInitialMesh()

  def preprocessTimeIter(self, timeIter):
    """
    Des.
    """

    self.SU2.PreprocessExtIter(timeIter)

  def remeshing(self):
    """
    Desctiption.
    """      
 
  def exit(self):
    """
    Description.
    """

    self.SU2.Postprocessing()

  def fakeSolidSolver(self, time):
    """
    Des.
    """


