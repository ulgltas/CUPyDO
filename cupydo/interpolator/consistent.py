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

interpolator.py
Interface interpolator and communicator.
Authors : David THOMAS, Marco Lucio CERQUAGLIA, Romain BOMAN

'''

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import numpy as np
import sys
import ccupydo
from ..utilities import *
from ..interfaceData import FlexInterfaceData
from ..interfaceData import InterfaceMatrix
from ..linearSolver import LinearSolver
from .interpolator import InterfaceInterpolator

np.set_printoptions(threshold=sys.maxsize)

# ----------------------------------------------------------------------
#    Consistent Interpolator class
# ----------------------------------------------------------------------

class ConsistentInterpolator(InterfaceInterpolator):
    def __init__(self, Manager, FluidSolver, SolidSolver, p, mpiComm = None):

        InterfaceInterpolator.__init__(self, Manager, FluidSolver, SolidSolver, p, mpiComm)
        mpiPrint('\nSetting non-matching consistent interpolator...', mpiComm)

        self.d = self.nDim+1
        self.SolverA = None
        self.SolverC = None

    def getLinearSolvers(self):

        return [self.SolverA, self.SolverC]

    def checkConservation(self):

        mpiPrint('No conservation check for consistent interpolation.', self.mpiComm)

    def generateInterfaceData(self):

        if self.manager.mechanical:
            self.prevSolidInterfaceDisplacement = FlexInterfaceData(self.ns + self.d, 3, self.mpiComm)
            self.solidInterfaceDisplacement = FlexInterfaceData(self.ns + self.d, 3, self.mpiComm)
            self.fluidInterfaceDisplacement = FlexInterfaceData(self.nf, 3, self.mpiComm)
            self.solidInterfaceLoads = FlexInterfaceData(self.ns, 6, self.mpiComm)
            self.fluidInterfaceLoads = FlexInterfaceData(self.nf + self.d, 6, self.mpiComm)
            if self.manager.computation == 'adjoint':
                self.solidInterfaceAdjointDisplacement = FlexInterfaceData(self.ns + self.d, 3, self.mpiComm)
                self.fluidInterfaceAdjointDisplacement = FlexInterfaceData(self.nf, 3, self.mpiComm)
                self.solidInterfaceAdjointLoads = FlexInterfaceData(self.ns, 6, self.mpiComm)
                self.fluidInterfaceAdjointLoads = FlexInterfaceData(self.nf + self.d, 6, self.mpiComm)

        if self.manager.thermal :
            if self.chtTransferMethod == 'TFFB':
                self.solidInterfaceTemperature = FlexInterfaceData(self.ns, 1, self.mpiComm)
                self.fluidInterfaceTemperature = FlexInterfaceData(self.nf + self.d, 1, self.mpiComm)
                self.solidInterfaceHeatFlux = FlexInterfaceData(self.ns + self.d, 3, self.mpiComm)
                self.fluidInterfaceHeatFlux = FlexInterfaceData(self.nf, 3, self.mpiComm)
            elif self.chtTransferMethod == 'FFTB':
                self.solidInterfaceTemperature = FlexInterfaceData(self.ns + self.d, 1, self.mpiComm)
                self.fluidInterfaceTemperature = FlexInterfaceData(self.nf, 1, self.mpiComm)
                self.solidInterfaceHeatFlux = FlexInterfaceData(self.ns, 3, self.mpiComm)
                self.fluidInterfaceHeatFlux = FlexInterfaceData(self.nf + self.d, 3, self.mpiComm)
                self.fluidInterfaceNormalHeatFlux = FlexInterfaceData(self.nf + self.d, 1, self.mpiComm)
                self.solidInterfaceNormalHeatFlux = FlexInterfaceData(self.ns, 1, self.mpiComm)
            elif self.chtTransferMethod == 'hFTB':
                self.fluidInterfaceRobinTemperature = FlexInterfaceData(self.nf + self.d, 1, self.mpiComm)
                self.solidInterfaceRobinTemperature = FlexInterfaceData(self.ns, 1, self.mpiComm)
                self.solidInterfaceTemperature = FlexInterfaceData(self.ns + self.d, 1, self.mpiComm)
                self.fluidInterfaceTemperature = FlexInterfaceData(self.nf, 1, self.mpiComm)
            elif self.chtTransferMethod == 'hFFB':
                self.fluidInterfaceRobinTemperature = FlexInterfaceData(self.nf + self.d, 1, self.mpiComm)
                self.solidInterfaceRobinTemperature = FlexInterfaceData(self.ns, 1, self.mpiComm)
                self.solidInterfaceHeatFlux = FlexInterfaceData(self.ns + self.d, 3, self.mpiComm)
                self.fluidInterfaceHeatFlux = FlexInterfaceData(self.nf, 3, self.mpiComm)

        self.A = InterfaceMatrix((self.ns+self.d,self.ns+self.d), self.mpiComm)
        self.B = InterfaceMatrix((self.nf,self.ns+self.d), self.mpiComm)
        self.C = InterfaceMatrix((self.nf+self.d,self.nf+self.d), self.mpiComm)
        self.D = InterfaceMatrix((self.ns,self.nf+self.d), self.mpiComm)

    def generateMapping(self):

        solidInterfaceProcessors = self.manager.getSolidInterfaceProcessors()
        fluidInterfaceProcessors = self.manager.getFluidInterfaceProcessors()
        solidPhysicalInterfaceNodesDistribution = self.manager.getSolidPhysicalInterfaceNodesDistribution()
        fluidPhysicalInterfaceNodesDistribution = self.manager.getFluidPhysicalInterfaceNodesDistribution()

        mpiPrint('\nBuilding interpolation matrices...', self.mpiComm)
        mpiPrint('\nBuilding matrix A of size {} X {}...'.format(self.ns, self.ns), self.mpiComm)
        
        # Fill the matrix A
        if self.mpiComm != None:
            for iProc in solidInterfaceProcessors:
                if self.myid == iProc:
                    localSolidInterface_array_X, localSolidInterface_array_Y, localSolidInterface_array_Z = self.SolidSolver.getNodalInitialPositions()
                    for jProc in solidInterfaceProcessors:
                        self.mpiComm.Isend(localSolidInterface_array_X, dest=jProc, tag=1)
                        self.mpiComm.Isend(localSolidInterface_array_Y, dest=jProc, tag=2)
                        self.mpiComm.Isend(localSolidInterface_array_Z, dest=jProc, tag=3)
                if self.myid in solidInterfaceProcessors:
                    sizeOfBuff = solidPhysicalInterfaceNodesDistribution[iProc]
                    solidInterfaceBuffRcv_X = np.zeros(sizeOfBuff)
                    solidInterfaceBuffRcv_Y = np.zeros(sizeOfBuff)
                    solidInterfaceBuffRcv_Z = np.zeros(sizeOfBuff)
                    req = self.mpiComm.Irecv(solidInterfaceBuffRcv_X, iProc, tag=1)
                    req.Wait()
                    req = self.mpiComm.Irecv(solidInterfaceBuffRcv_Y, iProc, tag=2)
                    req.Wait()
                    req = self.mpiComm.Irecv(solidInterfaceBuffRcv_Z, iProc, tag=3)
                    req.Wait()
                    self.fillMatrixA(solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc)
        else:
            localSolidInterface_array_X, localSolidInterface_array_Y, localSolidInterface_array_Z = self.SolidSolver.getNodalInitialPositions()
            self.fillMatrixA(localSolidInterface_array_X, localSolidInterface_array_Y, localSolidInterface_array_Z, 0)

        mpiBarrier(self.mpiComm)
        mpiPrint("\nAssembling A...", self.mpiComm)
        start = tm.time()
        self.A.assemble()
        mpiBarrier(self.mpiComm)
        stop = tm.time()
        mpiPrint('Assembly performed in {} s'.format(stop-start), self.mpiComm)
        mpiPrint('Matrix A is built.', self.mpiComm)
        mpiPrint('\nBuilding matrix B & D of size {} X {} & {} X {}...'.format(self.nf, self.ns, self.ns, self.nf), self.mpiComm)

        # Fill the matrix B & D
        if self.mpiComm != None:
            for iProc in solidInterfaceProcessors:
                if self.myid == iProc:
                    for jProc in fluidInterfaceProcessors:
                        self.mpiComm.Isend(localSolidInterface_array_X, dest=jProc, tag=1)
                        self.mpiComm.Isend(localSolidInterface_array_Y, dest=jProc, tag=2)
                        self.mpiComm.Isend(localSolidInterface_array_Z, dest=jProc, tag=3)
                if self.myid in fluidInterfaceProcessors:
                    sizeOfBuff = solidPhysicalInterfaceNodesDistribution[iProc]
                    solidInterfaceBuffRcv_X = np.zeros(sizeOfBuff)
                    solidInterfaceBuffRcv_Y = np.zeros(sizeOfBuff)
                    solidInterfaceBuffRcv_Z = np.zeros(sizeOfBuff)
                    req = self.mpiComm.Irecv(solidInterfaceBuffRcv_X, iProc, tag=1)
                    req.Wait()
                    req = self.mpiComm.Irecv(solidInterfaceBuffRcv_Y, iProc, tag=2)
                    req.Wait()
                    req = self.mpiComm.Irecv(solidInterfaceBuffRcv_Z, iProc, tag=3)
                    req.Wait()
                    self.fillMatrixBD(solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc)
        else:
            self.fillMatrixBD(localSolidInterface_array_X, localSolidInterface_array_Y, localSolidInterface_array_Z, 0)

        mpiBarrier(self.mpiComm)
        mpiPrint("\nAssembling B & D...", self.mpiComm)
        start = tm.time()
        self.B.assemble()
        mpiBarrier(self.mpiComm)
        self.D.assemble()
        mpiBarrier(self.mpiComm)
        stop = tm.time()
        mpiPrint('Assembly performed in {} s'.format(stop-start), self.mpiComm)
        mpiPrint('Matrix B & D are built.', self.mpiComm)
        mpiPrint('\nBuilding matrix C of size {} X {}...'.format(self.nf, self.nf), self.mpiComm)

        # Fill the matrix C
        if self.mpiComm != None:
            for iProc in fluidInterfaceProcessors:
                if self.myid == iProc:
                    localFluidInterface_array_X, localFluidInterface_array_Y, localFluidInterface_array_Z = self.FluidSolver.getNodalInitialPositions()
                    for jProc in fluidInterfaceProcessors:
                        self.mpiComm.Send(localFluidInterface_array_X, dest=jProc, tag=1)
                        self.mpiComm.Send(localFluidInterface_array_Y, dest=jProc, tag=2)
                        self.mpiComm.Send(localFluidInterface_array_Z, dest=jProc, tag=3)
                if self.myid in fluidInterfaceProcessors:
                    sizeOfBuff = fluidPhysicalInterfaceNodesDistribution[iProc]
                    fluidInterfaceBuffRcv_X = np.zeros(sizeOfBuff)
                    fluidInterfaceBuffRcv_Y = np.zeros(sizeOfBuff)
                    fluidInterfaceBuffRcv_Z = np.zeros(sizeOfBuff)
                    self.mpiComm.Recv(fluidInterfaceBuffRcv_X, iProc, tag=1)
                    self.mpiComm.Recv(fluidInterfaceBuffRcv_Y, iProc, tag=2)
                    self.mpiComm.Recv(fluidInterfaceBuffRcv_Z, iProc, tag=3)
                    self.fillMatrixC(fluidInterfaceBuffRcv_X, fluidInterfaceBuffRcv_Y, fluidInterfaceBuffRcv_Z, iProc)
        else:
            localFluidInterface_array_X, localFluidInterface_array_Y, localFluidInterface_array_Z = self.FluidSolver.getNodalInitialPositions()
            self.fillMatrixC(localFluidInterface_array_X, localFluidInterface_array_Y, localFluidInterface_array_Z, 0)

        mpiBarrier(self.mpiComm)
        mpiPrint("\nAssembling C...", self.mpiComm)
        start = tm.time()
        self.C.assemble()
        mpiBarrier(self.mpiComm)
        stop = tm.time()
        mpiPrint('Assembly performed in {} s'.format(stop-start), self.mpiComm)
        mpiPrint('Matrix C is built.', self.mpiComm)

        self.SolverA = LinearSolver(self.A, self.mpiComm)
        self.SolverC = LinearSolver(self.C, self.mpiComm)

    def interpolateFluidToSolid(self, fluidInterfaceData, solidInterfaceData):

        dim = fluidInterfaceData.getDim()
        gamma_array = FlexInterfaceData(self.nf + self.d, dim, self.mpiComm)

        self.SolverC.solve(fluidInterfaceData, gamma_array)
        self.D.mult(gamma_array, solidInterfaceData)

    def interpolateSolidToFluid(self, solidInterfaceData, fluidInterfaceData):

        dim = solidInterfaceData.getDim()
        gamma_array = FlexInterfaceData(self.ns + self.d, dim, self.mpiComm)

        self.SolverA.solve(solidInterfaceData, gamma_array)
        self.B.mult(gamma_array, fluidInterfaceData)

# ----------------------------------------------------------------------
#    RBF Consistent Interpolator class
# ----------------------------------------------------------------------

class ConsistentRBFInterpolator(ConsistentInterpolator):
    def __init__(self, Manager, FluidSolver, SolidSolver, p, mpiComm= None):

        ConsistentInterpolator.__init__(self, Manager, FluidSolver, SolidSolver, p, mpiComm)
        mpiPrint('\nSetting interpolation with Radial Basis Functions...', mpiComm)

        self.radius = p['rbfRadius']
        self.generateInterfaceData()
        self.generateMapping()

    def generateInterfaceData(self):

        ConsistentInterpolator.generateInterfaceData(self)
        mpiPrint('Generating interface data for consistent RBF interpolator...', self.mpiComm)

        self.A.createSparseFullAlloc()
        self.B.createSparseFullAlloc()
        self.C.createSparseFullAlloc()
        self.D.createSparseFullAlloc()

    def fillMatrixA(self, solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc):

        localSolidInterface_array_X_init, localSolidInterface_array_Y_init, localSolidInterface_array_Z_init = self.SolidSolver.getNodalInitialPositions()
        start = tm.time()
        ccupydo.CInterpolator.consistent_RBF_fillMatrixA(self, localSolidInterface_array_X_init, localSolidInterface_array_Y_init, localSolidInterface_array_Z_init,
                                              solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, self.A, iProc, 1.01*self.radius)
        stop = tm.time()
        print('Built A on rank {} in {} s'.format(self.myid,stop-start))

    def fillMatrixBD(self, solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc):

        localFluidInterface_array_X_init, localFluidInterface_array_Y_init, localFluidInterface_array_Z_init = self.FluidSolver.getNodalInitialPositions()
        start = tm.time()
        ccupydo.CInterpolator.consistent_RBF_fillMatrixBD(self, localFluidInterface_array_X_init, localFluidInterface_array_Y_init, localFluidInterface_array_Z_init,
                                              solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, self.B, self.D, iProc, 1.01*self.radius)
        stop = tm.time()
        print('Built B & D on rank {} in {} s'.format(self.myid,stop-start))

    def fillMatrixC(self, fluidInterfaceBuffRcv_X, fluidInterfaceBuffRcv_Y, fluidInterfaceBuffRcv_Z, iProc):

        localFluidInterface_array_X_init, localFluidInterface_array_Y_init, localFluidInterface_array_Z_init = self.FluidSolver.getNodalInitialPositions()
        start = tm.time()
        ccupydo.CInterpolator.consistent_RBF_fillMatrixC(self, localFluidInterface_array_X_init, localFluidInterface_array_Y_init, localFluidInterface_array_Z_init,
                                              fluidInterfaceBuffRcv_X, fluidInterfaceBuffRcv_Y, fluidInterfaceBuffRcv_Z, self.C, iProc, 1.01*self.radius)
        stop = tm.time()
        print('Built C on rank {} in {} s'.format(self.myid,stop-start))

# ----------------------------------------------------------------------
#    TPS Consistent Interpolator class
# ----------------------------------------------------------------------

class ConsistentTPSInterpolator(ConsistentInterpolator):
    def __init__(self, Manager, FluidSolver, SolidSolver, p, mpiComm= None):

        ConsistentInterpolator.__init__(self, Manager, FluidSolver, SolidSolver, p, mpiComm)
        mpiPrint('\nSetting consistent interpolation with Thin Plate Spline...', self.mpiComm)

        self.generateInterfaceData()
        self.generateMapping()

    def generateInterfaceData(self):

        ConsistentInterpolator.generateInterfaceData(self)
        mpiPrint('Generating interface data for consistent TPS interpolator...', self.mpiComm)

        self.A.createDense()
        self.B.createDense()
        self.C.createDense()
        self.D.createDense()

    def fillMatrixA(self, solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc):

        localSolidInterface_array_X_init, localSolidInterface_array_Y_init, localSolidInterface_array_Z_init = self.SolidSolver.getNodalInitialPositions()
        start = tm.time()
        ccupydo.CInterpolator.consistent_TPS_fillMatrixA(self, localSolidInterface_array_X_init, localSolidInterface_array_Y_init, localSolidInterface_array_Z_init,
                                              solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, self.A, iProc)
        stop = tm.time()
        print('Built A on rank {} in {} s'.format(self.myid,stop-start))

    def fillMatrixBD(self, solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc):

        localFluidInterface_array_X_init, localFluidInterface_array_Y_init, localFluidInterface_array_Z_init = self.FluidSolver.getNodalInitialPositions()
        start = tm.time()
        ccupydo.CInterpolator.consistent_TPS_fillMatrixBD(self, localFluidInterface_array_X_init, localFluidInterface_array_Y_init, localFluidInterface_array_Z_init,
                                              solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, self.B, self.D, iProc)
        stop = tm.time()
        print('Built B & D on rank {} in {} s'.format(self.myid,stop-start))

    def fillMatrixC(self, fluidInterfaceBuffRcv_X, fluidInterfaceBuffRcv_Y, fluidInterfaceBuffRcv_Z, iProc):

        localFluidInterface_array_X_init, localFluidInterface_array_Y_init, localFluidInterface_array_Z_init = self.FluidSolver.getNodalInitialPositions()
        start = tm.time()
        ccupydo.CInterpolator.consistent_TPS_fillMatrixC(self, localFluidInterface_array_X_init, localFluidInterface_array_Y_init, localFluidInterface_array_Z_init,
                                              fluidInterfaceBuffRcv_X, fluidInterfaceBuffRcv_Y, fluidInterfaceBuffRcv_Z, self.C, iProc)
        stop = tm.time()
        print('Built C on rank {} in {} s'.format(self.myid,stop-start))
