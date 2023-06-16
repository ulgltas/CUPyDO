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

np.set_printoptions(threshold=sys.maxsize)

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

        mpiPrint("Initializing FSI Interpolator",mpiComm,titlePrint)

        ccupydo.CInterpolator.__init__(self, Manager)

        self.manager = Manager
        self.SolidSolver = SolidSolver
        self.FluidSolver = FluidSolver

        self.mappingTimer = Timer()

        self.nf = int(self.manager.getNumberOfFluidInterfaceNodes())
        self.ns = int(self.manager.getNumberOfSolidInterfaceNodes())
        self.nf_loc = int(self.manager.getNumberOfLocalFluidInterfaceNodes())
        self.ns_loc = int(self.manager.getNumberOfLocalSolidInterfaceNodes())
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

        self.solidInterfaceAdjointDisplacement = None
        self.fluidInterfaceAdjointDisplacement = None
        self.solidInterfaceAdjointLoads = None
        self.fluidInterfaceAdjointLoads = None

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

    def getAdjointLoadsFromSolidSolver(self):
        """
        Des.
        """

        if self.myid in self.manager.getSolidInterfaceProcessors():
            localSolidInterfaceAdjointLoad_X, localSolidInterfaceAdjointLoad_Y, localSolidInterfaceAdjointLoad_Z = self.SolidSolver.getNodalAdjointLoads()
            for iVertex in range(self.ns_loc):
                iGlobalVertex = self.manager.getGlobalIndex('solid', self.myid, iVertex)
                self.solidInterfaceAdjointLoads[iGlobalVertex] = [localSolidInterfaceAdjointLoad_X[iVertex], localSolidInterfaceAdjointLoad_Y[iVertex], localSolidInterfaceAdjointLoad_Z[iVertex]]

        self.solidInterfaceAdjointLoads.assemble()

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
            localFluidInterfaceRobinTemperature = localFluidInterfaceTemperature - localFluidInterfaceNormalHeatFlux/self.heatTransferCoeff
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

    def getAdjointDisplacementFromFluidSolver(self):
        """
        Des.
        """

        if self.myid in self.manager.getFluidInterfaceProcessors():
            localFluidInterfaceAdjDisp_X, localFluidInterfaceAdjDisp_Y, localFluidInterfaceAdjDisp_Z = self.FluidSolver.getNodalAdjointDisplacement()
            for iVertex in range(self.nf_loc):
                iGlobalVertex = self.manager.getGlobalIndex('fluid', self.myid, iVertex)
                self.fluidInterfaceAdjointDisplacement[iGlobalVertex] = [localFluidInterfaceAdjDisp_X[iVertex], localFluidInterfaceAdjDisp_Y[iVertex], localFluidInterfaceAdjDisp_Z[iVertex]]

        self.fluidInterfaceAdjointDisplacement.assemble()

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
                    for key in list(fluidHaloNodesList[iProc].keys()):
                        globalIndex = fluidIndexing[key]
                        sendBuffHalo[key] = []
                        for iDim in range(fluidInterfaceData.nDim):
                            sendBuffHalo[key].append(fluidInterfaceData_array_recon[iDim][globalIndex])
                    iTagSend = 1
                    if iProc == 0: # In the master node use non-blocking and immediately receive them
                        for iDim in range(fluidInterfaceData.nDim):
                            self.mpiComm.Isend(sendBuff[iDim], dest=iProc, tag = iTagSend)
                            iTagSend += 1
                        sendBuffHalo_key = np.array(list(sendBuffHalo.keys()))
                        sendBuffHalo_values = np.empty((sendBuffHalo_key.size, 3),dtype=float)
                        for ii in range(sendBuffHalo_key.size):
                            sendBuffHalo_values[ii] = np.array(sendBuffHalo[sendBuffHalo_key[ii]])
                        self.mpiComm.Isend(np.array(sendBuffHalo_key.size), dest=iProc, tag=101)
                        self.mpiComm.Isend(sendBuffHalo_key, dest=iProc, tag=102)
                        self.mpiComm.Isend(sendBuffHalo_values, dest=iProc, tag=103)
                        localFluidInterfaceData_array = []
                        iTagRec = 1
                        for iDim in range(fluidInterfaceData.nDim):
                            local_array = np.zeros(self.nf_loc)
                            self.mpiComm.Recv(local_array, source=0, tag=iTagRec)
                            localFluidInterfaceData_array.append(local_array)
                            iTagRec += 1
                        nHaloNodesRcv = np.empty(1, dtype=int)
                        req = self.mpiComm.Irecv(nHaloNodesRcv, source=0, tag=101)
                        req.Wait()
                        rcvBuffHalo_keyBuff = np.empty(nHaloNodesRcv[0], dtype=int)
                        req = self.mpiComm.Irecv(rcvBuffHalo_keyBuff, source=0, tag=102)
                        req.Wait()
                        rcvBuffHalo_values = np.empty((nHaloNodesRcv[0],3), dtype=float)
                        req = self.mpiComm.Irecv(rcvBuffHalo_values, source=0, tag=103)
                        req.Wait()
                        for ii in range(len(rcvBuffHalo_keyBuff)):
                            haloNodesData_bis[rcvBuffHalo_keyBuff[ii]] = list(rcvBuffHalo_values[ii])
                        haloNodesData = haloNodesData_bis
                    else: # In other processors it's ok to send the buffers with blocking comms
                        for iDim in range(fluidInterfaceData.nDim):
                            self.mpiComm.Send(sendBuff[iDim], dest=iProc, tag = iTagSend)
                            iTagSend += 1
                        #self.mpiComm.send(sendBuffHalo, dest = iProc, tag=iTagSend)
                        sendBuffHalo_key = np.array(list(sendBuffHalo.keys()))
                        sendBuffHalo_values = np.empty((sendBuffHalo_key.size, 3),dtype=float)
                        for ii in range(sendBuffHalo_key.size):
                            sendBuffHalo_values[ii] = np.array(sendBuffHalo[sendBuffHalo_key[ii]])
                        self.mpiComm.Send(np.array(sendBuffHalo_key.size), dest=iProc, tag=101)
                        self.mpiComm.Send(sendBuffHalo_key, dest=iProc, tag=102)
                        self.mpiComm.Send(sendBuffHalo_values, dest=iProc, tag=103)
            elif self.myid in self.manager.getFluidInterfaceProcessors():
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
                localSolidInterfaceData_array = []
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
                    for key in list(solidHaloNodesList[iProc].keys()):
                        globalIndex = solidIndexing[key]
                        sendBuffHalo[key] = []
                        for iDim in range(solidInterfaceData.nDim):
                            sendBuffHalo[key].append(solidInterfaceData_array_recon[iDim][globalIndex])
                    iTagSend = 1
                    for iDim in range(solidInterfaceData.nDim):
                        self.mpiComm.Isend(sendBuff[iDim], dest=iProc, tag = iTagSend)
                        if iProc == 0:
                            local_array = np.zeros(self.ns_loc)
                            self.mpiComm.Recv(local_array, source=0, tag = iTagSend)
                            localSolidInterfaceData_array.append(local_array)
                        iTagSend += 1
                    #self.mpiComm.send(sendBuffHalo, dest = iProc, tag=iTagSend)
                    sendBuffHalo_key = np.array(list(sendBuffHalo.keys()))
                    sendBuffHalo_values = np.empty((sendBuffHalo_key.size, 3),dtype=float)
                    for ii in range(sendBuffHalo_key.size):
                        sendBuffHalo_values[ii] = np.array(sendBuffHalo[sendBuffHalo_key[ii]])
                    if iProc == 0:
                        self.mpiComm.Isend(np.array(sendBuffHalo_key.size), dest=iProc, tag=101)
                        self.mpiComm.Isend(sendBuffHalo_key, dest=iProc, tag=102)
                        self.mpiComm.Isend(sendBuffHalo_values, dest=iProc, tag=103)
                        nHaloNodesRcv = np.empty(1, dtype=int)
                        self.mpiComm.Recv(nHaloNodesRcv, source=0, tag=101)
                        rcvBuffHalo_keyBuff = np.empty(nHaloNodesRcv[0], dtype=int)
                        self.mpiComm.Recv(rcvBuffHalo_keyBuff, source=0, tag=102)
                        rcvBuffHalo_values = np.empty((nHaloNodesRcv[0],3), dtype=float)
                        self.mpiComm.Recv(rcvBuffHalo_values, source=0, tag=103)
                        for ii in range(len(rcvBuffHalo_keyBuff)):
                            haloNodesData_bis[rcvBuffHalo_keyBuff[ii]] = list(rcvBuffHalo_values[ii])
                        haloNodesData = haloNodesData_bis
                    else:
                        self.mpiComm.Send(np.array(sendBuffHalo_key.size), dest=iProc, tag=101)
                        self.mpiComm.Send(sendBuffHalo_key, dest=iProc, tag=102)
                        self.mpiComm.Send(sendBuffHalo_values, dest=iProc, tag=103)
            elif self.myid in self.manager.getSolidInterfaceProcessors():
                localSolidInterfaceData_array = []
                iTagRec = 1
                for iDim in range(solidInterfaceData.nDim):
                    local_array = np.zeros(self.ns_loc)
                    req = self.mpiComm.Irecv(local_array, source=0, tag = iTagRec)
                    req.Wait()
                    localSolidInterfaceData_array.append(local_array)
                    iTagRec += 1
                #haloNodesData = self.mpiComm.recv(source=0, tag=iTagRec)
                nHaloNodesRcv = np.empty(1, dtype=int)
                req = self.mpiComm.Irecv(nHaloNodesRcv, source=0, tag=101)
                req.Wait()
                rcvBuffHalo_keyBuff = np.empty(nHaloNodesRcv[0], dtype=int)
                req = self.mpiComm.Irecv(rcvBuffHalo_keyBuff, source=0, tag=102)
                req.Wait()
                rcvBuffHalo_values = np.empty((nHaloNodesRcv[0],3), dtype=float)
                req = self.mpiComm.Irecv(rcvBuffHalo_values, source=0, tag=103)
                req.Wait()
                for ii in range(len(rcvBuffHalo_keyBuff)):
                    haloNodesData_bis[rcvBuffHalo_keyBuff[ii]] = list(rcvBuffHalo_values[ii])
                haloNodesData = haloNodesData_bis

        return (localSolidInterfaceData_array, haloNodesData)

    def setLoadsToSolidSolver(self, dt):
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
                self.SolidSolver.applyNodalLoads(localSolidLoads_array[0], localSolidLoads_array[1], localSolidLoads_array[2], dt, haloNodesSolidLoads)
                FX = localSolidLoads_array[0].sum()
                FY = localSolidLoads_array[1].sum()
                FZ = localSolidLoads_array[2].sum()
            FXT = mpiAllReduce(self.mpiComm, FX)
            FYT = mpiAllReduce(self.mpiComm, FY)
            FZT = mpiAllReduce(self.mpiComm, FZ)
        else:
            self.SolidSolver.applyNodalLoads(self.solidInterfaceLoads.getDataArray(0), self.solidInterfaceLoads.getDataArray(1), self.solidInterfaceLoads.getDataArray(2), dt)
            FXT, FYT, FZT = self.solidInterfaceLoads.sum()

        mpiPrint("Checking f/s interface total force...", self.mpiComm)
        mpiPrint('Solid side (Fx, Fy, Fz) = ({}, {}, {})'.format(FXT, FYT, FZT), self.mpiComm)
        mpiPrint('Fluid side (Fx, Fy, Fz) = ({}, {}, {})'.format(FFX, FFY, FFZ), self.mpiComm)

    def setAdjointDisplacementToSolidSolver(self, dt):
        """
        Des.
        """

        # self.checkConservation()

        if self.mpiComm != None:
            (localSolidInterfaceAdjointDisplacement, haloNodesDisplacements) = self.redistributeDataToSolidSolver(self.solidInterfaceAdjointDisplacement)
#            mpiBarrier(self.mpiComm)
            if self.myid in self.manager.getSolidInterfaceProcessors():
                self.SolidSolver.applyNodalAdjointDisplacement(localSolidInterfaceAdjointDisplacement[0], localSolidInterfaceAdjointDisplacement[1], localSolidInterfaceAdjointDisplacement[2], haloNodesDisplacements, dt)
        else:
            self.SolidSolver.applyNodalAdjointDisplacement(self.solidInterfaceAdjointDisplacement.getDataArray(0), self.solidInterfaceAdjointDisplacement.getDataArray(1), self.solidInterfaceAdjointDisplacement.getDataArray(2), {}, dt)


    def setDisplacementToFluidSolver(self, dt):
        """
        Des.
        """

        self.checkConservation()

        if self.mpiComm != None:
            (localFluidInterfaceDisplacement, haloNodesDisplacements) = self.redistributeDataToFluidSolver(self.fluidInterfaceDisplacement)
#            mpiBarrier(self.mpiComm)
            if self.myid in self.manager.getFluidInterfaceProcessors():
                self.FluidSolver.applyNodalDisplacements(localFluidInterfaceDisplacement[0], localFluidInterfaceDisplacement[1], localFluidInterfaceDisplacement[2], localFluidInterfaceDisplacement[0], localFluidInterfaceDisplacement[1], localFluidInterfaceDisplacement[2], haloNodesDisplacements, dt)
        else:
            self.FluidSolver.applyNodalDisplacements(self.fluidInterfaceDisplacement.getDataArray(0), self.fluidInterfaceDisplacement.getDataArray(1), self.fluidInterfaceDisplacement.getDataArray(2), self.fluidInterfaceDisplacement.getDataArray(0), self.fluidInterfaceDisplacement.getDataArray(1), self.fluidInterfaceDisplacement.getDataArray(2), {}, dt)

    def setHeatFluxToFluidSolver(self, dt):
        """
        Description.
        """

        if self.mpiComm != None:
            (localFluidInterfaceHeatFlux, haloNodesHeatFlux) = self.redistributeDataToFluidSolver(self.fluidInterfaceHeatFlux)
            if self.myid in self.manager.getFluidInterfaceProcessors():
                self.FluidSolver.applyNodalHeatFluxes(localFluidInterfaceHeatFlux[0], localFluidInterfaceHeatFlux[1], localFluidInterfaceHeatFlux[2], dt)
        else:
            self.FluidSolver.applyNodalHeatFluxes(self.fluidInterfaceHeatFlux.getDataArray(0), self.fluidInterfaceHeatFlux.getDataArray(1), self.fluidInterfaceHeatFlux.getDataArray(2), dt)

    def setTemperatureToFluidSolver(self, dt):
        """
        Des.
        """

        if self.mpiComm != None:
            (localFluidInterfaceTemperature, haloNodesTemperature) = self.redistributeDataToFluidSolver(self.fluidInterfaceTemperature)
            if self.myid in self.manager.getFluidInterfaceProcessors():
                self.FluidSolver.applyNodalTemperatures(localFluidInterfaceTemperature[0], dt)
        else:
            self.FluidSolver.applyNodalTemperatures(self.fluidInterfaceTemperature.getDataArray(0), dt)

    def setAdjointLoadsToFluidSolver(self, dt):
        """
        Des.
        """

        # self.checkConservation()

        if self.mpiComm != None:
            (localFluidInterfaceAdjointLoad, haloNodesAdjointLoads) = self.redistributeDataToFluidSolver(self.fluidInterfaceAdjointLoads)
#            mpiBarrier(self.mpiComm)
            if self.myid in self.manager.getFluidInterfaceProcessors():
                self.FluidSolver.applyNodalAdjointLoads(localFluidInterfaceAdjointLoad[0], localFluidInterfaceAdjointLoad[1], localFluidInterfaceAdjointLoad[2], haloNodesAdjointLoads, dt)
        else:
            self.FluidSolver.applyNodalAdjointLoads(self.fluidInterfaceAdjointLoads.getDataArray(0), self.fluidInterfaceAdjointLoads.getDataArray(1), self.fluidInterfaceAdjointLoads.getDataArray(2), {}, dt)


    def setTemperatureToSolidSolver(self, dt):
        """
        Description
        """

        if self.mpiComm != None:
            (localSolidInterfaceTemperature, haloNodesTemperature) = self.redistributeDataToSolidSolver(self.solidInterfaceTemperature)
            if self.myid in self.manager.getSolidInterfaceProcessors():
                self.SolidSolver.applyNodalTemperatures(localSolidInterfaceTemperature[0], dt)
        else:
            self.SolidSolver.applyNodalTemperatures(self.solidInterfaceTemperature.getDataArray(0), dt)

    def setRobinHeatFluxToSolidSolver(self, dt):
        """
        Def
        """

        if self.mpiComm != None:
            (localSolidInterfaceRobinTemperature, haloNodesRobinTemperature) = self.redistributeDataToSolidSolver(self.solidInterfaceRobinTemperature)
            if self.myid in self.manager.getSolidInterfaceProcessors():
                localSolidInterfaceTemperature = self.SolidSolver.getNodalTemperatures()
                localSolidInterfaceRobinHeatFlux = self.heatTransferCoeff*(localSolidInterfaceTemperature-localSolidInterfaceRobinTemperature[0])
                self.SolidSolver.applyNodalNormalHeatFluxes(localSolidInterfaceRobinHeatFlux, dt)
        else:
            localSolidInterfaceTemperature = self.SolidSolver.getNodalTemperatures()
            localSolidInterfaceRobinHeatFlux = self.heatTransferCoeff*(localSolidInterfaceTemperature-self.solidInterfaceRobinTemperature.getDataArray(0), dt)
            self.SolidSolver.applyNodalNormalHeatFluxes(localSolidInterfaceRobinHeatFlux, dt)

    def setHeatFluxToSolidSolver(self, dt):
        """
        Des.
        """

        if self.mpiComm != None:
            (localSolidInterfaceNormalHeatFlux, haloNodesNormalHeatFlux) =  self.redistributeDataToSolidSolver(self.solidInterfaceNormalHeatFlux)
            if self.myid in self.manager.getSolidInterfaceProcessors():
                self.SolidSolver.applyNodalNormalHeatFluxes(localSolidInterfaceNormalHeatFlux[0], dt)
        else:
            self.SolidSolver.applyNodalNormalHeatFluxes(self.solidInterfaceNormalHeatFlux.getDataArray(0), dt)

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

    def interpolateSolidAdjointLoadsOnFluidMesh(self):
        """
        Description.
        """

        self.interpolateSolidToFluid(self.solidInterfaceAdjointLoads, self.fluidInterfaceAdjointLoads)

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

    def interpolateFluidAdjointDisplacementOnSolidMesh(self):
        """
        Des.
        """

        self.interpolateFluidToSolid(self.fluidInterfaceAdjointDisplacement, self.solidInterfaceAdjointDisplacement)


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
    