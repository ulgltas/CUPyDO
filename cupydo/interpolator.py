#! /usr/bin/env python
# -*- coding: latin-1; -*-

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

import ccupydo
from utilities import *
from interfaceData import FlexInterfaceData
from interfaceData import InterfaceMatrix
from linearSolver import LinearSolver

np.set_printoptions(threshold=np.nan)

# ----------------------------------------------------------------------
#    Interpolator class
# ----------------------------------------------------------------------

class InterfaceInterpolator(ccupydo.CInterpolator):
    """
    Interpolator of CUPyDO.
    Perform inteporlation of fluid-structure meshes.
    Inherited public members :
        -matching_fillMatrix()
        -TPS_fillMatrixA()
        -TPS_fillMatrixB()
        -RBF_fillMatrixA()
        -RBF_fillMatrixB()
        -PHI_TPS()
        -PHI_RBF()
        -distance()
    """

    def __init__(self, Manager, FluidSolver, SolidSolver, mpiComm = None, chtTransferMethod=None, heatTransferCoeff=1.0):
        """
        Description.
        """

        mpiPrint('\n***************************** Initializing FSI interpolator *****************************', mpiComm)

        ccupydo.CInterpolator.__init__(self, Manager)

        self.manager = Manager
        self.SolidSolver = SolidSolver
        self.FluidSolver = FluidSolver

        self.mappingTimer = Timer()

        self.nf = self.manager.getNumberOfFluidInterfaceNodes()
        self.ns = self.manager.getNumberOfSolidInterfaceNodes()
        self.nf_loc = self.manager.getNumberOfLocalFluidInterfaceNodes()
        self.ns_loc = self.manager.getNumberOfLocalSolidInterfaceNodes()
        self.nDim = self.manager.getnDim()

        self.d = 0

        if self.manager.thermal:
            self.chtTransferMethod = chtTransferMethod
            if self.chtTransferMethod not in ['TFFB','FFTB','hFTB','hFFB']:
                mpiPrint('CHT transfer method not specified or not recognized, using default TFFB',mpiComm)
                self.chtTransferMethod = 'TFFB'
        else:
            self.chtTransferMethod = None

        if self.chtTransferMethod in ['hFTB','hFFB']:
            self.heatTransferCoeff = heatTransferCoeff
        else:
            self.heatTransferCoeff = None

        self.mpiComm = mpiComm

        if self.mpiComm != None:
            self.myid = self.mpiComm.Get_rank()
            self.mpiSize = self.mpiComm.Get_size()
        else:
            self.myid = 0
            self.mpiSize = 1

        self.solidInterfaceDisplacement = None
        self.fluidInterfaceDisplacement = None
        self.solidInterfaceLoads = None
        self.fluidInterfaceLoads = None

        self.solidInterfaceHeatFlux = None
        self.fluidInterfaceHeatFlux = None
        self.solidInterfaceTemperature = None
        self.fluidInterfaceTemperature = None
        self.fluidInterfaceNormalHeatFlux = None
        self.solidInterfaceNormalHeatFlux = None
        self.fluidInterfaceRobinTemperature = None
        self.solidInterfaceRobinTemperature = None

    def checkTotalLoad(self):
        """
        Des.
        """

        FX, FY, FZ = self.solidInterfaceLoads.sum()

        FFX, FFY, FFZ = self.fluidInterfaceLoads.sum()

        mpiPrint("Checking f/s interface total force...", self.mpiComm)
        mpiPrint('Solid side (Fx, Fy, Fz) = ({}, {}, {})'.format(FX, FY, FZ), self.mpiComm)
        mpiPrint('Fluid side (Fx, Fy, Fz) = ({}, {}, {})'.format(FFX, FFY, FFZ), self.mpiComm)

    def getDisplacementFromSolidSolver(self):
        """
        Des.
        """

        if self.myid in self.manager.getSolidInterfaceProcessors():
            localSolidInterfaceDisp_X, localSolidInterfaceDisp_Y, localSolidInterfaceDisp_Z = self.SolidSolver.getNodalDisplacements()
            for iVertex in range(self.ns_loc):
                iGlobalVertex = self.manager.getGlobalIndex('solid', self.myid, iVertex)
                self.solidInterfaceDisplacement[iGlobalVertex] = [localSolidInterfaceDisp_X[iVertex], localSolidInterfaceDisp_Y[iVertex], localSolidInterfaceDisp_Z[iVertex]]

        self.solidInterfaceDisplacement.assemble()

    def getHeatFluxFromSolidSolver(self):
        """
        Des.
        """

        if self.myid in self.manager.getSolidInterfaceProcessors():
            localSolidInterfaceHeatFlux_X, localSolidInterfaceHeatFlux_Y, localSolidInterfaceHeatFlux_Z = self.SolidSolver.getNodalHeatFluxes()
            for iVertex in range(self.ns_loc):
                iGlobalVertex = self.manager.getGlobalIndex('solid', self.myid, iVertex)
                self.solidInterfaceHeatFlux[iGlobalVertex] = [localSolidInterfaceHeatFlux_X[iVertex], localSolidInterfaceHeatFlux_Y[iVertex], localSolidInterfaceHeatFlux_Z[iVertex]]

        self.solidInterfaceHeatFlux.assemble()

    def getLoadsFromFluidSolver(self):
        """
        Des.
        """

        if self.myid in self.manager.getFluidInterfaceProcessors():
            localFluidInterfaceLoad_X, localFluidInterfaceLoad_Y, localFluidInterfaceLoad_Z = self.FluidSolver.getNodalLoads()
            for iVertex in range(self.nf_loc):
                iGlobalVertex = self.manager.getGlobalIndex('fluid', self.myid, iVertex)
                self.fluidInterfaceLoads[iGlobalVertex] = [localFluidInterfaceLoad_X[iVertex], localFluidInterfaceLoad_Y[iVertex], localFluidInterfaceLoad_Z[iVertex]]

        self.fluidInterfaceLoads.assemble()

    def getTemperatureFromFluidSolver(self):
        """
        Des.
        """

        if self.myid in self.manager.getFluidInterfaceProcessors():
            localFluidInterfaceTemperature = self.FluidSolver.getNodalTemperatures()
            for iVertex in range(self.nf_loc):
                iGlobalVertex = self.manager.getGlobalIndex('fluid', self.myid, iVertex)
                self.fluidInterfaceTemperature[iGlobalVertex] = [localFluidInterfaceTemperature[iVertex]]

        self.fluidInterfaceTemperature.assemble()

    def getRobinTemperatureFromFluidSolver(self):
        """
        Des.
        """

        if self.myid in self.manager.getFluidInterfaceProcessors():
            localFluidInterfaceNormalHeatFlux = self.FluidSolver.getNodalNormalHeatFlux()
            localFluidInterfaceTemperature = self.FluidSolver.getNodalTemperatures()
            localFluidInterfaceRobinTemperature = localFluidInterfaceTemperature - (localFluidInterfaceNormalHeatFlux/self.heatTransferCoeff)
            for iVertex in range(self.nf_loc):
                iGlobalVertex = self.manager.getGlobalIndex('fluid', self.myid, iVertex)
                self.fluidInterfaceRobinTemperature[iGlobalVertex] = [localFluidInterfaceRobinTemperature[iVertex]]

        self.fluidInterfaceRobinTemperature.assemble()

    def getHeatFluxFromFluidSolver(self):
        """
        Des.
        """

        if self.myid in self.manager.getFluidInterfaceProcessors():
            localFluidInterfaceHeatFlux_X, localFluidInterfaceHeatFlux_Y, localFluidInterfaceHeatFlux_Z = self.FluidSolver.getNodalHeatFluxes()
            localFluidInterfaceNormalHeatFlux = self.FluidSolver.getNodalNormalHeatFlux()
            for iVertex in range(self.nf_loc):
                iGlobalVertex = self.manager.getGlobalIndex('fluid', self.myid, iVertex)
                self.fluidInterfaceHeatFlux[iGlobalVertex] = [localFluidInterfaceHeatFlux_X[iVertex], localFluidInterfaceHeatFlux_Y[iVertex], localFluidInterfaceHeatFlux_Z[iVertex]]
                self.fluidInterfaceNormalHeatFlux[iGlobalVertex] = [localFluidInterfaceNormalHeatFlux[iVertex]]

        self.fluidInterfaceHeatFlux.assemble()
        self.fluidInterfaceNormalHeatFlux.assemble()

    def redistributeDataToFluidSolver(self, fluidInterfaceData):
        """
        Description
        """

        localFluidInterfaceData_array = None
        haloNodesData = {}

        if self.mpiComm != None:
            localSize = fluidInterfaceData.getDataArray(0).shape[0]
            fluidInterfaceData_array_recon = []
            for iDim in range(fluidInterfaceData.nDim):
                array_recon = mpiGatherv(fluidInterfaceData.getDataArray(iDim), localSize, self.nf, self.mpiComm, 0)
                fluidInterfaceData_array_recon.append(array_recon)
            haloNodesData = {}
            haloNodesData_bis = {}
            if self.myid == 0:
                for iProc in self.manager.getFluidInterfaceProcessors():
                    fluidPhysicalInterfaceNodesDistribution = self.manager.getFluidPhysicalInterfaceNodesDistribution()
                    fluidGlobalIndexRange = self.manager.getFluidGlobalIndexRange()
                    sendBuff = []
                    for iDim in range(fluidInterfaceData.nDim):
                        sendBuff_i = np.zeros(fluidPhysicalInterfaceNodesDistribution[iProc])
                        sendBuff.append(sendBuff_i)
                    globalIndex = fluidGlobalIndexRange[iProc][0]
                    sendBuffHalo = {}
                    for iVertex in range(fluidPhysicalInterfaceNodesDistribution[iProc]):
                        for iDim in range(fluidInterfaceData.nDim):
                            sendBuff[iDim][iVertex] = fluidInterfaceData_array_recon[iDim][globalIndex]
                        globalIndex += 1
                    fluidHaloNodesList = self.manager.getFluidHaloNodesList()
                    fluidIndexing = self.manager.getFluidIndexing()
                    for key in fluidHaloNodesList[iProc].keys():
                        globalIndex = fluidIndexing[key]
                        sendBuffHalo[key] = []
                        for iDim in range(fluidInterfaceData.nDim):
                            sendBuffHalo[key].append(fluidInterfaceData_array_recon[iDim][globalIndex])
                    iTagSend = 1
                    for iDim in range(fluidInterfaceData.nDim):
                        self.mpiComm.Send(sendBuff[iDim], dest=iProc, tag = iTagSend)
                        iTagSend += 1
                    #self.mpiComm.send(sendBuffHalo, dest = iProc, tag=iTagSend)
                    sendBuffHalo_key = np.array(sendBuffHalo.keys())
                    sendBuffHalo_values = np.empty((sendBuffHalo_key.size, 3),dtype=float)
                    for ii in range(sendBuffHalo_key.size):
                        sendBuffHalo_values[ii] = np.array(sendBuffHalo[sendBuffHalo_key[ii]])
                    self.mpiComm.Send(np.array(sendBuffHalo_key.size), dest=iProc, tag=101)
                    self.mpiComm.Send(sendBuffHalo_key, dest=iProc, tag=102)
                    self.mpiComm.Send(sendBuffHalo_values, dest=iProc, tag=103)
            if self.myid in self.manager.getFluidInterfaceProcessors():
                localFluidInterfaceData_array = []
                iTagRec = 1
                for iDim in range(fluidInterfaceData.nDim):
                    local_array = np.zeros(self.nf_loc)
                    self.mpiComm.Recv(local_array, source=0, tag=iTagRec)
                    localFluidInterfaceData_array.append(local_array)
                    iTagRec += 1
                #haloNodesData = self.mpiComm.recv(source=0, tag=iTagRec)
                nHaloNodesRcv = np.empty(1, dtype=int)
                self.mpiComm.Recv(nHaloNodesRcv, source=0, tag=101)
                rcvBuffHalo_keyBuff = np.empty(nHaloNodesRcv[0], dtype=int)
                self.mpiComm.Recv(rcvBuffHalo_keyBuff, source=0, tag=102)
                rcvBuffHalo_values = np.empty((nHaloNodesRcv[0],3), dtype=float)
                self.mpiComm.Recv(rcvBuffHalo_values, source=0, tag=103)
                for ii in range(len(rcvBuffHalo_keyBuff)):
                    haloNodesData_bis[rcvBuffHalo_keyBuff[ii]] = list(rcvBuffHalo_values[ii])
                haloNodesData = haloNodesData_bis


        return (localFluidInterfaceData_array, haloNodesData)

    def redistributeDataToSolidSolver(self, solidInterfaceData):
        """
        Des.
        """

        localSolidInterfaceData_array = None
        haloNodesData = {}
        haloNodesData_bis = {}

        if self.mpiComm != None:
            localSize = solidInterfaceData.getDataArray(0).shape[0]
            solidInterfaceData_array_recon = []
            for iDim in range(solidInterfaceData.nDim):
                array_recon = mpiGatherv(solidInterfaceData.getDataArray(iDim), localSize, self.ns+self.d, self.mpiComm, 0)
                solidInterfaceData_array_recon.append(array_recon)
            haloNodesData = {}
            if self.myid == 0:
                for iProc in self.manager.getSolidInterfaceProcessors():
                    solidPhysicalInterfaceNodesDistribution = self.manager.getSolidPhysicalInterfaceNodesDistribution()
                    solidGlobalIndexRange = self.manager.getSolidGlobalIndexRange()
                    sendBuff = []
                    for iDim in range(solidInterfaceData.nDim):
                        sendBuff_i = np.zeros(solidPhysicalInterfaceNodesDistribution[iProc])
                        sendBuff.append(sendBuff_i)
                    globalIndex = solidGlobalIndexRange[iProc][0]
                    sendBuffHalo = {}
                    for iVertex in range(solidPhysicalInterfaceNodesDistribution[iProc]):
                        for iDim in range(solidInterfaceData.nDim):
                            sendBuff[iDim][iVertex] = solidInterfaceData_array_recon[iDim][globalIndex]
                        globalIndex += 1
                    solidHaloNodesList = self.manager.getSolidHaloNodesList()
                    solidIndexing = self.manager.getSolidIndexing()
                    for key in solidHaloNodesList[iProc].keys():
                        globalIndex = solidIndexing[key]
                        sendBuffHalo[key] = []
                        for iDim in range(solidInterfaceData.nDim):
                            sendBuffHalo[key].append(solidInterfaceData_array_recon[iDim][globalIndex])
                    iTagSend = 1
                    for iDim in range(solidInterfaceData.nDim):
                        self.mpiComm.Send(sendBuff[iDim], dest=iProc, tag = iTagSend)
                        iTagSend += 1
                    #self.mpiComm.send(sendBuffHalo, dest = iProc, tag=iTagSend)
                    sendBuffHalo_key = np.array(sendBuffHalo.keys())
                    sendBuffHalo_values = np.empty((sendBuffHalo_key.size, 3),dtype=float)
                    for ii in range(sendBuffHalo_key.size):
                        sendBuffHalo_values[ii] = np.array(sendBuffHalo[sendBuffHalo_key[ii]])
                    self.mpiComm.Send(np.array(sendBuffHalo_key.size), dest=iProc, tag=101)
                    self.mpiComm.Send(sendBuffHalo_key, dest=iProc, tag=102)
                    self.mpiComm.Send(sendBuffHalo_values, dest=iProc, tag=103)
            if self.myid in self.manager.getSolidInterfaceProcessors():
                localSolidInterfaceData_array = []
                iTagRec = 1
                for iDim in range(solidInterfaceData.nDim):
                    local_array = np.zeros(self.ns_loc)
                    self.mpiComm.Recv(local_array, source=0, tag = iTagRec)
                    localSolidInterfaceData_array.append(local_array)
                    iTagRec += 1
                #haloNodesData = self.mpiComm.recv(source=0, tag=iTagRec)
                nHaloNodesRcv = np.empty(1, dtype=int)
                self.mpiComm.Recv(nHaloNodesRcv, source=0, tag=101)
                rcvBuffHalo_keyBuff = np.empty(nHaloNodesRcv[0], dtype=int)
                self.mpiComm.Recv(rcvBuffHalo_keyBuff, source=0, tag=102)
                rcvBuffHalo_values = np.empty((nHaloNodesRcv[0],3), dtype=float)
                self.mpiComm.Recv(rcvBuffHalo_values, source=0, tag=103)
                for ii in range(len(rcvBuffHalo_keyBuff)):
                    haloNodesData_bis[rcvBuffHalo_keyBuff[ii]] = list(rcvBuffHalo_values[ii])
                haloNodesData = haloNodesData_bis

        return (localSolidInterfaceData_array, haloNodesData)

    def setLoadsToSolidSolver(self, time):
        """
        des.
        """

        FFX, FFY, FFZ = self.fluidInterfaceLoads.sum()


        FX = 0.
        FY = 0.
        FZ = 0.

        FXT = 0.
        FYT = 0.
        FZT = 0.

        if self.mpiComm != None:
            (localSolidLoads_array, haloNodesSolidLoads) = self.redistributeDataToSolidSolver(self.solidInterfaceLoads)
            if self.myid in self.manager.getSolidInterfaceProcessors():
                self.SolidSolver.applyNodalLoads(localSolidLoads_array[0], localSolidLoads_array[1], localSolidLoads_array[2], time)
                FX = localSolidLoads_array[0].sum()
                FY = localSolidLoads_array[1].sum()
                FZ = localSolidLoads_array[2].sum()
            FXT = mpiAllReduce(self.mpiComm, FX)
            FYT = mpiAllReduce(self.mpiComm, FY)
            FZT = mpiAllReduce(self.mpiComm, FZ)
        else:
            self.SolidSolver.applyNodalLoads(self.solidInterfaceLoads.getDataArray(0), self.solidInterfaceLoads.getDataArray(1), self.solidInterfaceLoads.getDataArray(2), time)
            FXT, FYT, FZT = self.solidInterfaceLoads.sum()

        mpiPrint("Checking f/s interface total force...", self.mpiComm)
        mpiPrint('Solid side (Fx, Fy, Fz) = ({}, {}, {})'.format(FXT, FYT, FZT), self.mpiComm)
        mpiPrint('Fluid side (Fx, Fy, Fz) = ({}, {}, {})'.format(FFX, FFY, FFZ), self.mpiComm)

    def setDisplacementToFluidSolver(self, time):
        """
        Des.
        """

        self.checkConservation()

        if self.mpiComm != None:
            (localFluidInterfaceDisplacement, haloNodesDisplacements) = self.redistributeDataToFluidSolver(self.fluidInterfaceDisplacement)
            if self.myid in self.manager.getFluidInterfaceProcessors():
                self.FluidSolver.applyNodalDisplacements(localFluidInterfaceDisplacement[0], localFluidInterfaceDisplacement[1], localFluidInterfaceDisplacement[2], localFluidInterfaceDisplacement[0], localFluidInterfaceDisplacement[1], localFluidInterfaceDisplacement[2], haloNodesDisplacements, time)
        else:
            self.FluidSolver.applyNodalDisplacements(self.fluidInterfaceDisplacement.getDataArray(0), self.fluidInterfaceDisplacement.getDataArray(1), self.fluidInterfaceDisplacement.getDataArray(2), self.fluidInterfaceDisplacement.getDataArray(0), self.fluidInterfaceDisplacement.getDataArray(1), self.fluidInterfaceDisplacement.getDataArray(2), {}, time)

    def setHeatFluxToFluidSolver(self, time):
        """
        Description.
        """

        if self.mpiComm != None:
            (localFluidInterfaceHeatFlux, haloNodesHeatFlux) = self.redistributeDataToFluidSolver(self.fluidInterfaceHeatFlux)
            if self.myid in self.manager.getFluidInterfaceProcessors():
                self.FluidSolver.applyNodalHeatFluxes(localFluidInterfaceHeatFlux[0], localFluidInterfaceHeatFlux[1], localFluidInterfaceHeatFlux[2], time)
        else:
            self.FluidSolver.applyNodalHeatFluxes(self.fluidInterfaceHeatFlux.getDataArray(0), self.fluidInterfaceHeatFlux.getDataArray(1), self.fluidInterfaceHeatFlux.getDataArray(2), time)

    def setTemperatureToFluidSolver(self, time):
        """
        Des.
        """

        if self.mpiComm != None:
            (localFluidInterfaceTemperature, haloNodesTemperature) = self.redistributeDataToFluidSolver(self.fluidInterfaceTemperature)
            if self.myid in self.manager.getFluidInterfaceProcessors():
                self.FluidSolver.applyNodalTemperatures(localFluidInterfaceTemperature[0], time)
        else:
            self.FluidSolver.applyNodalTemperatures(self.fluidInterfaceTemperature.getDataArray(0), time)

    def setTemperatureToSolidSolver(self, time):
        """
        Description
        """

        if self.mpiComm != None:
            (localSolidInterfaceTemperature, haloNodesTemperature) = self.redistributeDataToSolidSolver(self.solidInterfaceTemperature)
            if self.myid in self.manager.getSolidInterfaceProcessors():
                self.SolidSolver.applyNodalTemperatures(localSolidInterfaceTemperature[0], time)
        else:
            self.SolidSolver.applyNodalTemperatures(self.solidInterfaceTemperature.getDataArray(0), time)

    def setRobinHeatFluxToSolidSolver(self, time):
        """
        Def
        """

        if self.mpiComm != None:
            (localSolidInterfaceRobinTemperature, haloNodesRobinTemperature) = self.redistributeDataToSolidSolver(self.solidInterfaceRobinTemperature)
            if self.myid in self.manager.getSolidInterfaceProcessors():
                localSolidInterfaceTemperature = self.SolidSolver.getNodalTemperatures()
                localSolidInterfaceRobinHeatFlux = self.heatTransferCoeff*(localSolidInterfaceTemperature-localSolidInterfaceRobinTemperature[0])
                self.SolidSolver.applyNodalNormalHeatFluxes(localSolidInterfaceRobinHeatFlux, time)
        else:
            localSolidInterfaceTemperature = self.SolidSolver.getNodalTemperatures()
            localSolidInterfaceRobinHeatFlux = self.heatTransferCoeff*(localSolidInterfaceTemperature-self.solidInterfaceRobinTemperature.getDataArray(0), time)
            self.SolidSolver.applyNodalNormalHeatFluxes(localSolidInterfaceRobinHeatFlux, time)

    def setHeatFluxToSolidSolver(self, time):
        """
        Des.
        """

        if self.mpiComm != None:
            (localSolidInterfaceNormalHeatFlux, haloNodesNormalHeatFlux) =  self.redistributeDataToSolidSolver(self.solidInterfaceNormalHeatFlux)
            if self.myid in self.manager.getSolidInterfaceProcessors():
                self.SolidSolver.applyNodalNormalHeatFluxes(localSolidInterfaceNormalHeatFlux[0], time)
        else:
            self.SolidSolver.applyNodalNormalHeatFluxes(self.solidInterfaceNormalHeatFlux.getDataArray(0), time)

    def interpolateFluidLoadsOnSolidMesh(self):
        """
        Description
        """

        self.interpolateFluidToSolid(self.fluidInterfaceLoads, self.solidInterfaceLoads)

    def interpolateSolidDisplacementOnFluidMesh(self):
        """
        Description.
        """

        self.interpolateSolidToFluid(self.solidInterfaceDisplacement, self.fluidInterfaceDisplacement)

    def interpolateSolidHeatFluxOnFluidMesh(self):
        """
        Description.
        """

        self.interpolateSolidToFluid(self.solidInterfaceHeatFlux, self.fluidInterfaceHeatFlux)


    def interpolateSolidTemperatureOnFluidMesh(self):
        """
        Description
        """

        self.interpolateSolidToFluid(self.solidInterfaceTemperature, self.fluidInterfaceTemperature)

    def interpolateFluidHeatFluxOnSolidMesh(self):
        """
        Description.
        """

        self.interpolateFluidToSolid(self.fluidInterfaceHeatFlux, self.solidInterfaceHeatFlux)
        self.interpolateFluidToSolid(self.fluidInterfaceNormalHeatFlux, self.solidInterfaceNormalHeatFlux)

    def interpolateFluidTemperatureOnSolidMesh(self):
        """
        Description.
        """

        self.interpolateFluidToSolid(self.fluidInterfaceTemperature, self.solidInterfaceTemperature)

    def interpolateFluidRobinTemperatureOnSolidMesh(self):
        """
        Des.
        """

        self.interpolateFluidToSolid(self.fluidInterfaceRobinTemperature, self.solidInterfaceRobinTemperature)

    def getNs(self):
        """
        Des.
        """

        return self.ns

    def getNf(self):
        """
        Des.
        """

        return self.nf

    def getd(self):
        """
        Des.
        """

        return self.d

class MatchingMeshesInterpolator(InterfaceInterpolator):
    """
    Description.
    """

    def __init__(self, Manager, FluidSolver, SolidSolver, mpiComm = None, chtTransferMethod=None, heatTransferCoeff=1.0):
        """
        Description
        """

        InterfaceInterpolator.__init__(self, Manager, FluidSolver, SolidSolver, mpiComm, chtTransferMethod, heatTransferCoeff)

        mpiPrint('\nSetting matching meshes interpolator...', mpiComm)

        if self.nf != self.ns:
            raise Exception("Fluid and solid interface must have the same number of nodes for matching meshes ! ")
        ccupydo.CInterpolator.matching_initSearch(self)

        self.generateInterfaceData()

        self.generateMapping()

    def checkConservation(self):
        """
        Des.
        """

        WSX, WSY, WSZ = self.solidInterfaceLoads.dot(self.solidInterfaceDisplacement)

        WFX, WFY, WFZ = self.fluidInterfaceLoads.dot(self.fluidInterfaceDisplacement)

        mpiPrint("Checking f/s interface conservation...", self.mpiComm)
        mpiPrint('Solid side (Wx, Wy, Wz) = ({}, {}, {})'.format(WSX, WSY, WSZ), self.mpiComm)
        mpiPrint('Fluid side (Wx, Wy, Wz) = ({}, {}, {})'.format(WFX, WFY, WFZ), self.mpiComm)

    def generateInterfaceData(self):
        """
        Des.
        """

        if self.manager.mechanical:
            self.solidInterfaceDisplacement = FlexInterfaceData(self.ns, 3, self.mpiComm)
            self.fluidInterfaceDisplacement = FlexInterfaceData(self.nf, 3, self.mpiComm)
            self.solidInterfaceLoads = FlexInterfaceData(self.ns, 3, self.mpiComm)
            self.fluidInterfaceLoads = FlexInterfaceData(self.nf, 3, self.mpiComm)

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
        """
        Des.
        """

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
                        self.mpiComm.Send(localSolidInterface_array_X, dest=jProc, tag=1)
                        self.mpiComm.Send(localSolidInterface_array_Y, dest=jProc, tag=2)
                        self.mpiComm.Send(localSolidInterface_array_Z, dest=jProc, tag=3)
                if self.myid in fluidInterfaceProcessors:
                    sizeOfBuff = solidPhysicalInterfaceNodesDistribution[iProc]
                    solidInterfaceBuffRcv_X = np.zeros(sizeOfBuff)
                    solidInterfaceBuffRcv_Y = np.zeros(sizeOfBuff)
                    solidInterfaceBuffRcv_Z = np.zeros(sizeOfBuff)
                    self.mpiComm.Recv(solidInterfaceBuffRcv_X, iProc, tag=1)
                    self.mpiComm.Recv(solidInterfaceBuffRcv_Y, iProc, tag=2)
                    self.mpiComm.Recv(solidInterfaceBuffRcv_Z, iProc, tag=3)
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
        """
        Des.
        """

        localFluidInterface_array_X_init, localFluidInterface_array_Y_init, localFluidInterface_array_Z_init = self.FluidSolver.getNodalInitialPositions()

        print('Matching mapping search on rank {}...'.format(self.myid))
        start = tm.time()
        ccupydo.CInterpolator.matching_search(self, localFluidInterface_array_X_init, localFluidInterface_array_Y_init, localFluidInterface_array_Z_init,
                                              solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc)
        stop = tm.time()
        print('Search on rank {} in {} s'.format(self.myid,stop-start))

    def fillMatrix(self):
        """
        Des.
        """

        print('Building H on rank {}...'.format(self.myid))
        start = tm.time()
        ccupydo.CInterpolator.matching_fillMatrix(self, self.H, self.H_T)
        stop = tm.time()
        print('Built H on rank {} in {} s'.format(self.myid,stop-start))

    def interpolateFluidToSolid(self, fluidInterfaceData, solidInterfaceData):
        """
        des.
        """

        self.H_T.mult(fluidInterfaceData, solidInterfaceData)

    def interpolateSolidToFluid(self, solidInterfaceData, fluidInterfaceData):
        """
        Des.
        """

        self.H.mult(solidInterfaceData, fluidInterfaceData)

class ConservativeInterpolator(InterfaceInterpolator):
    """
    Description.
    """

    def __init__(self, Manager, FluidSolver, SolidSolver, mpiComm = None, chtTransferMethod=None, heatTransferCoeff=1.0):
        """
        Des.
        """

        InterfaceInterpolator.__init__(self, Manager, FluidSolver, SolidSolver, mpiComm, chtTransferMethod, heatTransferCoeff)

        mpiPrint('\nSetting non-matching conservative interpolator...', mpiComm)

        self.d = self.nDim+1
        self.SolverA = None
        self.SolverA_T = None

    def getLinearSolvers(self):
        """
        Des.
        """

        return [self.SolverA, self.SolverA_T]

    def checkConservation(self):
        """
        Des.
        """

        WSX, WSY, WSZ = self.solidInterfaceLoads.dot(self.solidInterfaceDisplacement)

        WFX, WFY, WFZ = self.fluidInterfaceLoads.dot(self.fluidInterfaceDisplacement)

        mpiPrint("Checking f/s interface conservation...", self.mpiComm)
        mpiPrint('Solid side (Wx, Wy, Wz) = ({}, {}, {})'.format(WSX, WSY, WSZ), self.mpiComm)
        mpiPrint('Fluid side (Wx, Wy, Wz) = ({}, {}, {})'.format(WFX, WFY, WFZ), self.mpiComm)

    def generateInterfaceData(self):
        """
        Description.
        """

        if self.manager.mechanical:
            self.solidInterfaceDisplacement = FlexInterfaceData(self.ns + self.d, 3, self.mpiComm)
            self.fluidInterfaceDisplacement = FlexInterfaceData(self.nf, 3, self.mpiComm)
            self.solidInterfaceLoads = FlexInterfaceData(self.ns + self.d, 3, self.mpiComm)
            self.fluidInterfaceLoads = FlexInterfaceData(self.nf, 3, self.mpiComm)

        if self.manager.thermal :
            if self.chtTransferMethod == 'TFFB':
                self.solidInterfaceTemperature = FlexInterfaceData(self.ns + self.d, 1, self.mpiComm)
                self.fluidInterfaceTemperature = FlexInterfaceData(self.nf, 1, self.mpiComm)
                self.solidInterfaceHeatFlux = FlexInterfaceData(self.ns + self.d, 3, self.mpiComm)
                self.fluidInterfaceHeatFlux = FlexInterfaceData(self.nf, 3, self.mpiComm)
            elif self.chtTransferMethod == 'FFTB':
                self.solidInterfaceTemperature = FlexInterfaceData(self.ns + self.d, 1, self.mpiComm)
                self.fluidInterfaceTemperature = FlexInterfaceData(self.nf, 1, self.mpiComm)
                self.solidInterfaceHeatFlux = FlexInterfaceData(self.ns + self.d, 3, self.mpiComm)
                self.fluidInterfaceHeatFlux = FlexInterfaceData(self.nf, 3, self.mpiComm)
                self.fluidInterfaceNormalHeatFlux = FlexInterfaceData(self.nf, 1, self.mpiComm)
                self.solidInterfaceNormalHeatFlux = FlexInterfaceData(self.ns, 1, self.mpiComm)
            elif self.chtTransferMethod == 'hFTB':
                self.fluidInterfaceRobinTemperature = FlexInterfaceData(self.nf, 1, self.mpiComm)
                self.solidInterfaceRobinTemperature = FlexInterfaceData(self.ns, 1, self.mpiComm)
                self.solidInterfaceTemperature = FlexInterfaceData(self.ns + self.d, 1, self.mpiComm)
                self.fluidInterfaceTemperature = FlexInterfaceData(self.nf, 1, self.mpiComm)
            elif self.chtTransferMethod == 'hFFB':
                self.fluidInterfaceRobinTemperature = FlexInterfaceData(self.nf, 1, self.mpiComm)
                self.solidInterfaceRobinTemperature = FlexInterfaceData(self.ns, 1, self.mpiComm)
                self.solidInterfaceHeatFlux = FlexInterfaceData(self.ns + self.d, 3, self.mpiComm)
                self.fluidInterfaceHeatFlux = FlexInterfaceData(self.nf, 3, self.mpiComm)

        self.A = InterfaceMatrix((self.ns+self.d,self.ns+self.d), self.mpiComm)
        self.A_T = InterfaceMatrix((self.ns+self.d,self.ns+self.d), self.mpiComm)
        self.B = InterfaceMatrix((self.nf,self.ns+self.d), self.mpiComm)
        self.B_T = InterfaceMatrix((self.ns+self.d,self.nf), self.mpiComm)

    def generateMapping(self):
        """
        Des.
        """

        solidInterfaceProcessors = self.manager.getSolidInterfaceProcessors()
        fluidInterfaceProcessors = self.manager.getFluidInterfaceProcessors()
        solidPhysicalInterfaceNodesDistribution = self.manager.getSolidPhysicalInterfaceNodesDistribution()

        mpiPrint('\nBuilding interpolation matrices...', self.mpiComm)

        mpiPrint('\nBuilding matrix A of size {} X {}...'.format(self.ns, self.ns), self.mpiComm)
        # Fill the matrix A
        if self.mpiComm != None:
            for iProc in solidInterfaceProcessors:
                if self.myid == iProc:
                    localSolidInterface_array_X, localSolidInterface_array_Y, localSolidInterface_array_Z = self.SolidSolver.getNodalInitialPositions()
                    for jProc in solidInterfaceProcessors:
                        self.mpiComm.Send(localSolidInterface_array_X, dest=jProc, tag=1)
                        self.mpiComm.Send(localSolidInterface_array_Y, dest=jProc, tag=2)
                        self.mpiComm.Send(localSolidInterface_array_Z, dest=jProc, tag=3)
                if self.myid in solidInterfaceProcessors:
                    sizeOfBuff = solidPhysicalInterfaceNodesDistribution[iProc]
                    solidInterfaceBuffRcv_X = np.zeros(sizeOfBuff)
                    solidInterfaceBuffRcv_Y = np.zeros(sizeOfBuff)
                    solidInterfaceBuffRcv_Z = np.zeros(sizeOfBuff)
                    self.mpiComm.Recv(solidInterfaceBuffRcv_X, iProc, tag=1)
                    self.mpiComm.Recv(solidInterfaceBuffRcv_Y, iProc, tag=2)
                    self.mpiComm.Recv(solidInterfaceBuffRcv_Z, iProc, tag=3)
                    self.fillMatrixA(solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc)
        else:
            localSolidInterface_array_X, localSolidInterface_array_Y, localSolidInterface_array_Z = self.SolidSolver.getNodalInitialPositions()
            self.fillMatrixA(localSolidInterface_array_X, localSolidInterface_array_Y, localSolidInterface_array_Z, 0)

        mpiBarrier(self.mpiComm)
        mpiPrint("\nAssembling A & A_T...", self.mpiComm)
        start = tm.time()
        self.A.assemble()
        mpiBarrier(self.mpiComm)
        self.A_T.assemble()
        mpiBarrier(self.mpiComm)
        stop = tm.time()
        mpiPrint('Assembly performed in {} s'.format(stop-start), self.mpiComm)
        mpiPrint('Matrix A is built.', self.mpiComm)

        mpiPrint('\nBuilding matrix B of size {} X {}...'.format(self.nf, self.ns), self.mpiComm)
        # Fill the matrix B
        if self.mpiComm != None:
            for iProc in solidInterfaceProcessors:
                if self.myid == iProc:
                    for jProc in fluidInterfaceProcessors:
                        self.mpiComm.Send(localSolidInterface_array_X, dest=jProc, tag=1)
                        self.mpiComm.Send(localSolidInterface_array_Y, dest=jProc, tag=2)
                        self.mpiComm.Send(localSolidInterface_array_Z, dest=jProc, tag=3)
                if self.myid in fluidInterfaceProcessors:
                    sizeOfBuff = solidPhysicalInterfaceNodesDistribution[iProc]
                    solidInterfaceBuffRcv_X = np.zeros(sizeOfBuff)
                    solidInterfaceBuffRcv_Y = np.zeros(sizeOfBuff)
                    solidInterfaceBuffRcv_Z = np.zeros(sizeOfBuff)
                    self.mpiComm.Recv(solidInterfaceBuffRcv_X, iProc, tag=1)
                    self.mpiComm.Recv(solidInterfaceBuffRcv_Y, iProc, tag=2)
                    self.mpiComm.Recv(solidInterfaceBuffRcv_Z, iProc, tag=3)
                    self.fillMatrixB(solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc)
        else:
            self.fillMatrixB(localSolidInterface_array_X, localSolidInterface_array_Y, localSolidInterface_array_Z, 0)

        mpiBarrier(self.mpiComm)
        mpiPrint("\nAssembling B & B_T...", self.mpiComm)
        start = tm.time()
        self.B.assemble()
        mpiBarrier(self.mpiComm)
        self.B_T.assemble()
        mpiBarrier(self.mpiComm)
        stop = tm.time()
        mpiPrint('Assembly performed in {} s'.format(stop-start), self.mpiComm)
        mpiPrint('Matrix B is built.', self.mpiComm)

        self.SolverA = LinearSolver(self.A, self.mpiComm)
        self.SolverA_T = LinearSolver(self.A_T, self.mpiComm)

    def interpolateFluidToSolid(self, fluidInterfaceData, solidInterfaceData):
        """
        des.
        """

        dim = fluidInterfaceData.getDim()
        gamma_array = FlexInterfaceData(self.ns + self.d, dim, self.mpiComm)

        self.B_T.mult(fluidInterfaceData, gamma_array)
        self.SolverA_T.solve(gamma_array, solidInterfaceData)

    def interpolateSolidToFluid(self, solidInterfaceData, fluidInterfaceData):
        """
        Des.
        """

        dim = solidInterfaceData.getDim()
        gamma_array = FlexInterfaceData(self.ns + self.d, dim, self.mpiComm)

        self.SolverA.solve(solidInterfaceData, gamma_array)
        self.B.mult(gamma_array, fluidInterfaceData)


class ConsistentInterpolator(InterfaceInterpolator):
    """
    Description.
    """

    def __init__(self, Manager, FluidSolver, SolidSolver, mpiComm = None, chtTransferMethod=None, heatTransferCoeff=1.0):
        """
        Des.
        """

        InterfaceInterpolator.__init__(self, Manager, FluidSolver, SolidSolver, mpiComm, chtTransferMethod, heatTransferCoeff)

        mpiPrint('\nSetting non-matching consistent interpolator...', mpiComm)

        self.d = self.nDim+1
        self.SolverA = None
        self.SolverC = None

    def getLinearSolvers(self):
        """
        Des.
        """

        return [self.SolverA, self.SolverC]

    def checkConservation(self):
        """
        Des.
        """

        mpiPrint('No conservation check for consistent interpolation.', self.mpiComm)

    def generateInterfaceData(self):
        """
        Description.
        """

        if self.manager.mechanical:
            self.solidInterfaceDisplacement = FlexInterfaceData(self.ns + self.d, 3, self.mpiComm)
            self.fluidInterfaceDisplacement = FlexInterfaceData(self.nf, 3, self.mpiComm)
            self.solidInterfaceLoads = FlexInterfaceData(self.ns, 3, self.mpiComm)
            self.fluidInterfaceLoads = FlexInterfaceData(self.nf + self.d, 3, self.mpiComm)

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
        """
        Des.
        """

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
                        self.mpiComm.Send(localSolidInterface_array_X, dest=jProc, tag=1)
                        self.mpiComm.Send(localSolidInterface_array_Y, dest=jProc, tag=2)
                        self.mpiComm.Send(localSolidInterface_array_Z, dest=jProc, tag=3)
                if self.myid in solidInterfaceProcessors:
                    sizeOfBuff = solidPhysicalInterfaceNodesDistribution[iProc]
                    solidInterfaceBuffRcv_X = np.zeros(sizeOfBuff)
                    solidInterfaceBuffRcv_Y = np.zeros(sizeOfBuff)
                    solidInterfaceBuffRcv_Z = np.zeros(sizeOfBuff)
                    self.mpiComm.Recv(solidInterfaceBuffRcv_X, iProc, tag=1)
                    self.mpiComm.Recv(solidInterfaceBuffRcv_Y, iProc, tag=2)
                    self.mpiComm.Recv(solidInterfaceBuffRcv_Z, iProc, tag=3)
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
                        self.mpiComm.Send(localSolidInterface_array_X, dest=jProc, tag=1)
                        self.mpiComm.Send(localSolidInterface_array_Y, dest=jProc, tag=2)
                        self.mpiComm.Send(localSolidInterface_array_Z, dest=jProc, tag=3)
                if self.myid in fluidInterfaceProcessors:
                    sizeOfBuff = solidPhysicalInterfaceNodesDistribution[iProc]
                    solidInterfaceBuffRcv_X = np.zeros(sizeOfBuff)
                    solidInterfaceBuffRcv_Y = np.zeros(sizeOfBuff)
                    solidInterfaceBuffRcv_Z = np.zeros(sizeOfBuff)
                    self.mpiComm.Recv(solidInterfaceBuffRcv_X, iProc, tag=1)
                    self.mpiComm.Recv(solidInterfaceBuffRcv_Y, iProc, tag=2)
                    self.mpiComm.Recv(solidInterfaceBuffRcv_Z, iProc, tag=3)
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
        """
        des.
        """

        dim = fluidInterfaceData.getDim()
        gamma_array = FlexInterfaceData(self.nf + self.d, dim, self.mpiComm)

        self.SolverC.solve(fluidInterfaceData, gamma_array)
        self.D.mult(gamma_array, solidInterfaceData)

    def interpolateSolidToFluid(self, solidInterfaceData, fluidInterfaceData):
        """
        Des.
        """

        dim = solidInterfaceData.getDim()
        gamma_array = FlexInterfaceData(self.ns + self.d, dim, self.mpiComm)

        self.SolverA.solve(solidInterfaceData, gamma_array)
        self.B.mult(gamma_array, fluidInterfaceData)

class RBFInterpolator(ConservativeInterpolator):
    """
    Description.
    """

    def __init__(self, Manager, FluidSolver, SolidSolver, RBFradius=0.1, mpiComm = None, chtTransferMethod=None, heatTransferCoeff=1.0):
        """"
        Description.
        """

        ConservativeInterpolator.__init__(self, Manager, FluidSolver, SolidSolver, mpiComm, chtTransferMethod, heatTransferCoeff)

        mpiPrint('\nSetting interpolation with Radial Basis Functions...', mpiComm)

        self.radius = RBFradius

        self.generateInterfaceData()

        self.generateMapping()


    def generateInterfaceData(self):
        """
        Des.
        """

        ConservativeInterpolator.generateInterfaceData(self)

        mpiPrint('Generating interface data for conservative RBF interpolator...', self.mpiComm)

        self.A.createSparseFullAlloc()
        self.A_T.createSparseFullAlloc()
        self.B.createSparseFullAlloc()
        self.B_T.createSparseFullAlloc()

    def fillMatrixA(self, solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc):
        """
        Description.
        """

        localSolidInterface_array_X_init, localSolidInterface_array_Y_init, localSolidInterface_array_Z_init = self.SolidSolver.getNodalInitialPositions()
        start = tm.time()
        ccupydo.CInterpolator.RBF_fillMatrixA(self, localSolidInterface_array_X_init, localSolidInterface_array_Y_init, localSolidInterface_array_Z_init,
                                              solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, self.A, self.A_T, iProc, 1.01*self.radius)
        stop = tm.time()
        print('Built A on rank {} in {} s'.format(self.myid,stop-start))


    def fillMatrixB(self, solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc):
        """
        Description.
        """

        localFluidInterface_array_X_init, localFluidInterface_array_Y_init, localFluidInterface_array_Z_init = self.FluidSolver.getNodalInitialPositions()
        start = tm.time()
        ccupydo.CInterpolator.RBF_fillMatrixB(self, localFluidInterface_array_X_init, localFluidInterface_array_Y_init, localFluidInterface_array_Z_init,
                                              solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, self.B, self.B_T, iProc, 1.01*self.radius)
        stop = tm.time()
        print('Built B on rank {} in {} s'.format(self.myid,stop-start))



class ConsistentRBFInterpolator(ConsistentInterpolator):
    """
    Description.
    """

    def __init__(self, Manager, FluidSolver, SolidSolver, RBFradius = 0.1, mpiComm= None, chtTransferMethod=None, heatTransferCoeff=1.0):
        """
        Des.
        """

        ConsistentInterpolator.__init__(self, Manager, FluidSolver, SolidSolver, mpiComm, chtTransferMethod, heatTransferCoeff)

        mpiPrint('\nSetting interpolation with Radial Basis Functions...', mpiComm)

        self.radius = RBFradius

        self.generateInterfaceData()

        self.generateMapping()

    def generateInterfaceData(self):
        """
        Des.
        """

        ConsistentInterpolator.generateInterfaceData(self)

        mpiPrint('Generating interface data for consistent RBF interpolator...', self.mpiComm)

        self.A.createSparseFullAlloc()
        self.B.createSparseFullAlloc()
        self.C.createSparseFullAlloc()
        self.D.createSparseFullAlloc()


    def fillMatrixA(self, solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc):
        """
        Description.
        """

        localSolidInterface_array_X_init, localSolidInterface_array_Y_init, localSolidInterface_array_Z_init = self.SolidSolver.getNodalInitialPositions()
        start = tm.time()
        ccupydo.CInterpolator.consistent_RBF_fillMatrixA(self, localSolidInterface_array_X_init, localSolidInterface_array_Y_init, localSolidInterface_array_Z_init,
                                              solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, self.A, iProc, 1.01*self.radius)
        stop = tm.time()
        print('Built A on rank {} in {} s'.format(self.myid,stop-start))

    def fillMatrixBD(self, solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc):
        """
        Description.
        """

        localFluidInterface_array_X_init, localFluidInterface_array_Y_init, localFluidInterface_array_Z_init = self.FluidSolver.getNodalInitialPositions()
        start = tm.time()
        ccupydo.CInterpolator.consistent_RBF_fillMatrixBD(self, localFluidInterface_array_X_init, localFluidInterface_array_Y_init, localFluidInterface_array_Z_init,
                                              solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, self.B, self.D, iProc, 1.01*self.radius)
        stop = tm.time()
        print('Built B & D on rank {} in {} s'.format(self.myid,stop-start))

    def fillMatrixC(self, fluidInterfaceBuffRcv_X, fluidInterfaceBuffRcv_Y, fluidInterfaceBuffRcv_Z, iProc):
        """
        Description.
        """

        localFluidInterface_array_X_init, localFluidInterface_array_Y_init, localFluidInterface_array_Z_init = self.FluidSolver.getNodalInitialPositions()
        start = tm.time()
        ccupydo.CInterpolator.consistent_RBF_fillMatrixC(self, localFluidInterface_array_X_init, localFluidInterface_array_Y_init, localFluidInterface_array_Z_init,
                                              fluidInterfaceBuffRcv_X, fluidInterfaceBuffRcv_Y, fluidInterfaceBuffRcv_Z, self.C, iProc, 1.01*self.radius)
        stop = tm.time()
        print('Built C on rank {} in {} s'.format(self.myid,stop-start))

class TPSInterpolator(ConservativeInterpolator):
    """
    Des.
    """

    def __init__(self, Manager, FluidSolver, SolidSolver, mpiComm=None, chtTransferMethod=None, heatTransferCoeff=1.0):
        """
        des.
        """

        ConservativeInterpolator.__init__(self, Manager, FluidSolver, SolidSolver, mpiComm, chtTransferMethod, heatTransferCoeff)

        mpiPrint('\nSetting interpolation with Thin Plate Spline...', self.mpiComm)

        self.generateInterfaceData()

        self.generateMapping()

    def generateInterfaceData(self):
        """
        Des.
        """

        ConservativeInterpolator.generateInterfaceData(self)

        mpiPrint('Generating interface data for TPS interpolator...', self.mpiComm)

        self.A.createDense()
        self.A_T.createDense()
        self.B.createDense()
        self.B_T.createDense()

    def fillMatrixA(self, solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc):
        """
        Description.
        """

        localSolidInterface_array_X_init, localSolidInterface_array_Y_init, localSolidInterface_array_Z_init = self.SolidSolver.getNodalInitialPositions()

        start = tm.time()
        ccupydo.CInterpolator.TPS_fillMatrixA(self, localSolidInterface_array_X_init, localSolidInterface_array_Y_init, localSolidInterface_array_Z_init,
                                              solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, self.A, self.A_T, iProc)
        stop = tm.time()
        print('Built A on rank {} in {} s'.format(self.myid,stop-start))

    def fillMatrixB(self, solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc):
        """
        Description.
        """

        localFluidInterface_array_X_init, localFluidInterface_array_Y_init, localFluidInterface_array_Z_init = self.FluidSolver.getNodalInitialPositions()

        start = tm.time()
        ccupydo.CInterpolator.TPS_fillMatrixB(self, localFluidInterface_array_X_init, localFluidInterface_array_Y_init, localFluidInterface_array_Z_init,
                                              solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, self.B, self.B_T, iProc)
        stop = tm.time()
        print('Built B on rank {} in {} s'.format(self.myid,stop-start))

class ConsistentTPSInterpolator(ConsistentInterpolator):
    """
    Description.
    """

    def __init__(self, Manager, FluidSolver, SolidSolver, mpiComm= None, chtTransferMethod=None, heatTransferCoeff=1.0):
        """
        Des.
        """

        ConsistentInterpolator.__init__(self, Manager, FluidSolver, SolidSolver, mpiComm, chtTransferMethod, heatTransferCoeff)

        mpiPrint('\nSetting consistent interpolation with Thin Plate Spline...', self.mpiComm)

        self.generateInterfaceData()

        self.generateMapping()

    def generateInterfaceData(self):
        """
        Des.
        """

        ConsistentInterpolator.generateInterfaceData(self)

        mpiPrint('Generating interface data for consistent TPS interpolator...', self.mpiComm)

        self.A.createDense()
        self.B.createDense()
        self.C.createDense()
        self.D.createDense()

    def fillMatrixA(self, solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc):
        """
        Des.
        """

        localSolidInterface_array_X_init, localSolidInterface_array_Y_init, localSolidInterface_array_Z_init = self.SolidSolver.getNodalInitialPositions()
        start = tm.time()
        ccupydo.CInterpolator.consistent_TPS_fillMatrixA(self, localSolidInterface_array_X_init, localSolidInterface_array_Y_init, localSolidInterface_array_Z_init,
                                              solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, self.A, iProc)
        stop = tm.time()
        print('Built A on rank {} in {} s'.format(self.myid,stop-start))

    def fillMatrixBD(self, solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc):
        """
        des.
        """

        localFluidInterface_array_X_init, localFluidInterface_array_Y_init, localFluidInterface_array_Z_init = self.FluidSolver.getNodalInitialPositions()
        start = tm.time()
        ccupydo.CInterpolator.consistent_TPS_fillMatrixBD(self, localFluidInterface_array_X_init, localFluidInterface_array_Y_init, localFluidInterface_array_Z_init,
                                              solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, self.B, self.D, iProc)
        stop = tm.time()
        print('Built B & D on rank {} in {} s'.format(self.myid,stop-start))

    def fillMatrixC(self, fluidInterfaceBuffRcv_X, fluidInterfaceBuffRcv_Y, fluidInterfaceBuffRcv_Z, iProc):
        """
        Des.
        """

        localFluidInterface_array_X_init, localFluidInterface_array_Y_init, localFluidInterface_array_Z_init = self.FluidSolver.getNodalInitialPositions()
        start = tm.time()
        ccupydo.CInterpolator.consistent_TPS_fillMatrixC(self, localFluidInterface_array_X_init, localFluidInterface_array_Y_init, localFluidInterface_array_Z_init,
                                              fluidInterfaceBuffRcv_X, fluidInterfaceBuffRcv_Y, fluidInterfaceBuffRcv_Z, self.C, iProc)
        stop = tm.time()
        print('Built C on rank {} in {} s'.format(self.myid,stop-start))
