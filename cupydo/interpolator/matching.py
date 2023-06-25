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
from .interpolator import InterfaceInterpolator

np.set_printoptions(threshold=sys.maxsize)

# ----------------------------------------------------------------------
#    Matching Mesh Interpolator class
# ----------------------------------------------------------------------

class MatchingMeshesInterpolator(InterfaceInterpolator):
    def __init__(self, Manager, FluidSolver, SolidSolver, mpiComm = None, chtTransferMethod=None, heatTransferCoeff=1.0):

        InterfaceInterpolator.__init__(self, Manager, FluidSolver, SolidSolver, mpiComm, chtTransferMethod, heatTransferCoeff)

        mpiPrint('\nSetting matching meshes interpolator...', mpiComm)

        if self.nf != self.ns:
            raise Exception("Fluid and solid interface must have the same number of nodes for matching meshes ! ")
        ccupydo.CInterpolator.matching_initSearch(self)

        self.generateInterfaceData()
        self.generateMapping()

    def checkConservation(self):

        WSX, WSY, WSZ = self.solidInterfaceLoads.dot(self.solidInterfaceDisplacement)
        WFX, WFY, WFZ = self.fluidInterfaceLoads.dot(self.fluidInterfaceDisplacement)

        mpiPrint("Checking f/s interface conservation...", self.mpiComm)
        mpiPrint('Solid side (Wx, Wy, Wz) = ({}, {}, {})'.format(WSX, WSY, WSZ), self.mpiComm)
        mpiPrint('Fluid side (Wx, Wy, Wz) = ({}, {}, {})'.format(WFX, WFY, WFZ), self.mpiComm)

    def generateInterfaceData(self):

        if self.manager.mechanical:
            self.prevSolidInterfaceDisplacement = FlexInterfaceData(self.ns, 3, self.mpiComm)
            self.solidInterfaceDisplacement = FlexInterfaceData(self.ns, 3, self.mpiComm)
            self.fluidInterfaceDisplacement = FlexInterfaceData(self.nf, 3, self.mpiComm)
            self.solidInterfaceLoads = FlexInterfaceData(self.ns, 3, self.mpiComm)
            self.fluidInterfaceLoads = FlexInterfaceData(self.nf, 3, self.mpiComm)
            if self.manager.computation == 'adjoint':
                self.solidInterfaceAdjointDisplacement = FlexInterfaceData(self.ns, 3, self.mpiComm)
                self.fluidInterfaceAdjointDisplacement = FlexInterfaceData(self.nf, 3, self.mpiComm)
                self.solidInterfaceAdjointLoads = FlexInterfaceData(self.ns, 3, self.mpiComm)
                self.fluidInterfaceAdjointLoads = FlexInterfaceData(self.nf, 3, self.mpiComm)

        if self.manager.thermal :
            if self.chtTransferMethod == 'TFFB':
                self.solidInterfaceTemperature = FlexInterfaceData(self.ns, 1, self.mpiComm)
                self.fluidInterfaceTemperature = FlexInterfaceData(self.nf, 1, self.mpiComm)
                self.solidInterfaceHeatFlux = FlexInterfaceData(self.ns, 3, self.mpiComm)
                self.fluidInterfaceHeatFlux = FlexInterfaceData(self.nf, 3, self.mpiComm)
            elif self.chtTransferMethod == 'FFTB':
                self.solidInterfaceTemperature = FlexInterfaceData(self.ns, 1, self.mpiComm)
                self.fluidInterfaceTemperature = FlexInterfaceData(self.nf, 1, self.mpiComm)
                self.solidInterfaceHeatFlux = FlexInterfaceData(self.ns, 3, self.mpiComm)
                self.fluidInterfaceHeatFlux = FlexInterfaceData(self.nf, 3, self.mpiComm)
                self.fluidInterfaceNormalHeatFlux = FlexInterfaceData(self.nf, 1, self.mpiComm)
                self.solidInterfaceNormalHeatFlux = FlexInterfaceData(self.ns, 1, self.mpiComm)
            elif self.chtTransferMethod == 'hFTB':
                self.fluidInterfaceRobinTemperature = FlexInterfaceData(self.nf, 1, self.mpiComm)
                self.solidInterfaceRobinTemperature = FlexInterfaceData(self.ns, 1, self.mpiComm)
                self.solidInterfaceTemperature = FlexInterfaceData(self.ns, 1, self.mpiComm)
                self.fluidInterfaceTemperature = FlexInterfaceData(self.nf, 1, self.mpiComm)
            elif self.chtTransferMethod == 'hFFB':
                self.fluidInterfaceRobinTemperature = FlexInterfaceData(self.nf, 1, self.mpiComm)
                self.solidInterfaceRobinTemperature = FlexInterfaceData(self.ns, 1, self.mpiComm)
                self.solidInterfaceHeatFlux = FlexInterfaceData(self.ns, 3, self.mpiComm)
                self.fluidInterfaceHeatFlux = FlexInterfaceData(self.nf, 3, self.mpiComm)

        self.H = InterfaceMatrix((self.nf,self.ns), self.mpiComm)
        self.H_T = InterfaceMatrix((self.ns,self.nf), self.mpiComm)
        self.H.createSparse(1,1)
        self.H_T.createSparse(1,1)

    def generateMapping(self):

        solidInterfaceProcessors = self.manager.getSolidInterfaceProcessors()
        fluidInterfaceProcessors = self.manager.getFluidInterfaceProcessors()
        solidPhysicalInterfaceNodesDistribution = self.manager.getSolidPhysicalInterfaceNodesDistribution()

        mpiPrint('\nBuilding interpolation matrix...', self.mpiComm)
        mpiPrint('\nBuilding matrix H of size {} X {}...'.format(self.nf, self.ns), self.mpiComm)
        self.mappingTimer.start()

        if self.mpiComm != None:
            for iProc in solidInterfaceProcessors:
                if self.myid == iProc:
                    localSolidInterface_array_X, localSolidInterface_array_Y, localSolidInterface_array_Z = self.SolidSolver.getNodalInitialPositions()
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
                    self.mappingSearch(solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc)
            if self.myid in fluidInterfaceProcessors:
                self.fillMatrix()
        else:
            localSolidInterface_array_X, localSolidInterface_array_Y, localSolidInterface_array_Z = self.SolidSolver.getNodalInitialPositions()
            self.mappingSearch(localSolidInterface_array_X, localSolidInterface_array_Y, localSolidInterface_array_Z, 0)
            self.fillMatrix()

        mpiBarrier(self.mpiComm)
        mpiPrint("\nAssembling H & H_T...", self.mpiComm)
        start = tm.time()
        self.H.assemble()
        mpiBarrier(self.mpiComm)
        self.H_T.assemble()
        mpiBarrier(self.mpiComm)
        stop = tm.time()
        mpiPrint('Assembly performed in {} s'.format(stop-start), self.mpiComm)
        mpiPrint('Matrix H is built.', self.mpiComm)

        self.mappingTimer.stop()
        self.mappingTimer.cumul()

    def mappingSearch(self, solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc):

        localFluidInterface_array_X_init, localFluidInterface_array_Y_init, localFluidInterface_array_Z_init = self.FluidSolver.getNodalInitialPositions()

        print('Matching mapping search on rank {}...'.format(self.myid))
        start = tm.time()
        ccupydo.CInterpolator.matching_search(self, localFluidInterface_array_X_init, localFluidInterface_array_Y_init, localFluidInterface_array_Z_init,
                                              solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc)
        stop = tm.time()
        print('Search on rank {} in {} s'.format(self.myid,stop-start))

    def fillMatrix(self):

        print('Building H on rank {}...'.format(self.myid))
        start = tm.time()
        ccupydo.CInterpolator.matching_fillMatrix(self, self.H, self.H_T)
        stop = tm.time()
        print('Built H on rank {} in {} s'.format(self.myid,stop-start))

    def interpolateFluidToSolid(self, fluidInterfaceData, solidInterfaceData):

        self.H_T.mult(fluidInterfaceData, solidInterfaceData)

    def interpolateSolidToFluid(self, solidInterfaceData, fluidInterfaceData):

        self.H.mult(solidInterfaceData, fluidInterfaceData)
