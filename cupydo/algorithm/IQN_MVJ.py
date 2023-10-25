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

algorithm.py
Defines the coupling algorithms of CUPyDO.
Authors : David THOMAS, Marco Lucio CERQUAGLIA, Romain BOMAN

'''

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

from math import *
import numpy as np
import sys
import sys

from ..utilities import *
from ..interfaceData import FlexInterfaceData
from .BGS import AlgorithmBGSStaticRelax

np.set_printoptions(threshold=sys.maxsize)

# ----------------------------------------------------------------------
#    IQN MVJ Algorithm class
# ----------------------------------------------------------------------

class AlgorithmIQN_MVJ(AlgorithmBGSStaticRelax):
    def __init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, Criterion, p, mpiComm):
        AlgorithmBGSStaticRelax.__init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, Criterion, p, mpiComm)

        self.resetInternalVars()

        # --- Tolerance and type of filtering : None, Degroote2, and Haelterman.
        self.tollQR = 1.0e-1
        self.qrFilter = 'Haelterman'

        # --- Option which allows to build the tangent matrix of a given time step using differences with respect to the first FSI iteration (delta_r_k = r_k+1 - r_0) instead of the previous iteration (delta_r_k = r_k+1 - r_k) --- #
        self.computeTangentMatrixBasedOnFirstIt = p['firstItTgtMat']

    def filter(self, V, W):

        if self.qrFilter == None: # No filtering is applied to V and W
            return V, W

        if self.qrFilter == 'Degroote2': # QR filtering as described by J. Degroote et al. CMAME, 199, 2085-2098 (2010).
            Q, R, V, W = QRfiltering(V, W, self.tollQR)
            return V, W
        
        elif self.qrFilter == 'Haelterman': # 'Modified' QR filtering as described by R. Haelterman et al. Computers and Structures, 171, 9-17 (2016).
            Q, R, V, W = QRfiltering_mod(V, W, self.tollQR)
            return V, W
        
        else:
            raise NameError('IQN-MVJ Algorithm: the QR filtering technique is unknown!')
        
    def resetInternalVars(self):

        ns = self.interfaceInterpolator.getNs()

        # --- Global V and W matrices for IQN-MVJ algorithm --- #
        if self.manager.mechanical:
            self.invJprev = np.zeros((self.manager.nDim*ns,self.manager.nDim*ns))
            self.invJ = np.zeros((self.manager.nDim*ns,self.manager.nDim*ns))
        if self.manager.thermal:
            self.invJprevCHT = np.zeros((self.manager.nDim*ns,self.manager.nDim*ns))
            self.invJCHT = np.zeros((self.manager.nDim*ns,self.manager.nDim*ns))

        # --- Indicate if a BGS iteration must be performed --- #
        if self.manager.mechanical:
            self.makeBGS = True
        if self.manager.thermal:
            self.makeBGS_CHT = True

    def fsiCoupling(self):
        """
        Interface Quasi Newton - Multi Vector Jacobian (IQN-MVJ) method for strong coupling FSI
        """

        mpiPrint("Enter IQN-MNJ Strong Coupling FSI",self.mpiComm,titlePrint)

        self.FSIIter = 0
        ns = self.interfaceInterpolator.getNs()
        d = self.interfaceInterpolator.getd()

        # --- Initialize all the quantities used in the IQN-MVJ method --- #
        if self.manager.mechanical:
            self.solidInterfaceResidual0 = FlexInterfaceData(ns+d, 3, self.mpiComm)
            self.solidInterfaceDisplacement_tilde = FlexInterfaceData(ns, 3, self.mpiComm)
            self.solidInterfaceDisplacement_tilde1 = FlexInterfaceData(ns, 3, self.mpiComm)
        if self.manager.thermal:
            self.solidInterfaceResidual0_CHT = FlexInterfaceData(ns+d, 1, self.mpiComm)
            self.solidInterfaceTemperature_tilde = FlexInterfaceData(ns, 1, self.mpiComm)
            self.solidInterfaceTemperature_tilde1 = FlexInterfaceData(ns, 1, self.mpiComm)

        # --- Initialises V, W --- #
        if self.manager.mechanical:
            self.Vk = []
            self.Wk = []
        if self.manager.thermal:
            self.VkCHT = []
            self.WkCHT = []

        while (self.FSIIter < self.nbFSIIterMax):
            mpiPrint("\n>>>> FSI iteration {} <<<<\n".format(self.FSIIter), self.mpiComm)

            # --- Solid to fluid mechanical transfer --- #
            if self.manager.mechanical:
                self.solidToFluidMechaTransfer()
                mpiPrint('\nPerforming mesh deformation...\n', self.mpiComm)
                self.meshDefTimer.start()
                self.FluidSolver.meshUpdate(self.step.timeIter)
                self.meshDefTimer.stop()
                self.meshDefTimer.cumul()

            # --- Solid to fluid thermal transfer --- #
            if self.manager.thermal:
                self.solidToFluidThermalTransfer()
                
            self.FluidSolver.boundaryConditionsUpdate()

            # --- Fluid solver call for FSI subiteration --- #
            mpiPrint('\nLaunching fluid solver...', self.mpiComm)
            self.fluidSolverTimer.start()
            verif = self.FluidSolver.run(*self.step.timeFrame())
            self.fluidSolverTimer.stop()
            self.fluidSolverTimer.cumul()
            mpiBarrier(self.mpiComm)

            # --- Check if the fluid solver succeeded --- #
            if not verif:
                self.resetInternalVars()
                return False

            # --- Fluid to solid mechanical transfer --- #
            if self.manager.mechanical:
                mpiPrint('\nProcessing interface fluid loads...\n', self.mpiComm)
                self.fluidToSolidMechaTransfer()

            if self.manager.thermal:
                # --- Fluid to solid thermal transfer --- #
                mpiPrint('\nProcessing interface thermal quantities...\n', self.mpiComm)
                self.fluidToSolidThermalTransfer()

            mpiBarrier(self.mpiComm)

            # --- Solid solver call for FSI subiteration --- #
            mpiPrint('\nLaunching solid solver...\n', self.mpiComm)
            if self.myid in self.manager.getSolidSolverProcessors():
                self.solidSolverTimer.start()
                verif = self.SolidSolver.run(*self.step.timeFrame())
                self.solidSolverTimer.stop()
                self.solidSolverTimer.cumul()

            # --- Check if the solid solver succeeded --- #
            try: solidProc = int(self.manager.getSolidSolverProcessors())
            except: raise Exception('Only one solid solver process is supported yet')
            verif = mpiScatter(verif, self.mpiComm, solidProc)
            self.solidHasRun = True

            if not verif:
                self.resetInternalVars()
                return False

            if self.manager.mechanical:
                # --- Compute and monitor the FSI residual --- #
                self.computeSolidInterfaceResidual()
                mpiPrint('\nFSI error value : {}\n'.format(self.criterion.epsilon), self.mpiComm)
                self.relaxMultiVectorJacobian()

            if self.manager.thermal:
                # --- Compute and monitor the FSI residual --- #
                self.computeSolidInterfaceResidual_CHT()
                mpiPrint('\nCHT error value : {}\n'.format(self.criterion.epsilonCHT), self.mpiComm)
                self.relaxMultiVectorJacobian_CHT()

            # --- Update the FSI iteration and history --- #
            if self.writeInFSIloop == True:
                self.writeRealTimeData()
            self.FSIIter += 1

            # --- Compute and monitor the FSI residual then update the Jacobian --- #
            if self.criterion.isVerified():
                mpiPrint("IQN-MVJ is Converged",self.mpiComm,titlePrint)
                if self.manager.mechanical:
                    self.invJprev = np.copy(self.invJ)
                if self.manager.thermal:
                    self.invJprevCHT = np.copy(self.invJCHT)
                return True
        
        # --- Reset the Jacobians because the coupling did not converge --- #
        self.resetInternalVars()
        return False


    def relaxMultiVectorJacobian(self):

        d = self.interfaceInterpolator.getd()
        ns = self.interfaceInterpolator.getNs()
        delta_ds = FlexInterfaceData(ns+d, 3, self.mpiComm)

        # --- Initialize d_tilde for the construction of the Wk matrix -- #
        if self.myid in self.manager.getSolidInterfaceProcessors():
            localSolidInterfaceDisp_X, localSolidInterfaceDisp_Y, localSolidInterfaceDisp_Z = self.SolidSolver.getNodalDisplacements()
            for iVertex in range(self.manager.getNumberOfLocalSolidInterfaceNodes()):
                iGlobalVertex = self.manager.getGlobalIndex('solid', self.myid, iVertex)
                self.solidInterfaceDisplacement_tilde[iGlobalVertex] = [localSolidInterfaceDisp_X[iVertex], localSolidInterfaceDisp_Y[iVertex], localSolidInterfaceDisp_Z[iVertex]]

        self.solidInterfaceDisplacement_tilde.assemble()
        
        # --- Relax the solid position --- #
        if self.makeBGS:
            self.invJprev = np.zeros((self.manager.nDim*ns,self.manager.nDim*ns))
            self.invJ = np.zeros((self.manager.nDim*ns,self.manager.nDim*ns))
            self.relaxSolidPosition()
            self.makeBGS = False

        # --- Construct Vk and Wk matrices for the computation of the approximated tangent matrix --- #
        else:
            mpiPrint('\nCorrect solid interface displacements using IQN-MVJ method...\n', self.mpiComm)
            
            # --- Start gathering on root process --- #
            res_X_Gat, res_Y_Gat, res_Z_Gat = mpiGatherInterfaceData(self.solidInterfaceResidual, ns+d, self.mpiComm, 0)
            solidInterfaceResidual0_X_Gat, solidInterfaceResidual0_Y_Gat, solidInterfaceResidual0_Z_Gat = mpiGatherInterfaceData(self.solidInterfaceResidual0, ns+d, self.mpiComm, 0)
            solidInterfaceDisplacement_tilde_X_Gat, solidInterfaceDisplacement_tilde_Y_Gat, solidInterfaceDisplacement_tilde_Z_Gat = mpiGatherInterfaceData(self.solidInterfaceDisplacement_tilde, ns, self.mpiComm, 0)
            solidInterfaceDisplacement_tilde1_X_Gat, solidInterfaceDisplacement_tilde1_Y_Gat, solidInterfaceDisplacement_tilde1_Z_Gat = mpiGatherInterfaceData(self.solidInterfaceDisplacement_tilde1, ns, self.mpiComm, 0)

            # --- Copies for operating on, length=ns, not ns+d --- #
            if self.myid == 0:
                res_X_Gat_C = res_X_Gat[:ns]
                res_Y_Gat_C = res_Y_Gat[:ns]
                res_Z_Gat_C = res_Z_Gat[:ns]
                solidInterfaceResidual0_X_Gat_C = solidInterfaceResidual0_X_Gat[:ns]
                solidInterfaceResidual0_Y_Gat_C = solidInterfaceResidual0_Y_Gat[:ns]
                solidInterfaceResidual0_Z_Gat_C = solidInterfaceResidual0_Z_Gat[:ns]


                # --- Use J from the previous time step because Vk and Wk are empty --- #
                if self.FSIIter == 0:
                    self.invJ = self.invJprev.copy()

                # -- Vk and Wk matrices are enriched only starting from the second iteration --- #
                else:
                    if self.manager.nDim == 3:
                        delta_res = np.concatenate([res_X_Gat_C - solidInterfaceResidual0_X_Gat_C, res_Y_Gat_C - solidInterfaceResidual0_Y_Gat_C, res_Z_Gat_C - solidInterfaceResidual0_Z_Gat_C], axis=0)
                        delta_d = np.concatenate([solidInterfaceDisplacement_tilde_X_Gat - solidInterfaceDisplacement_tilde1_X_Gat, solidInterfaceDisplacement_tilde_Y_Gat - solidInterfaceDisplacement_tilde1_Y_Gat, solidInterfaceDisplacement_tilde_Z_Gat - solidInterfaceDisplacement_tilde1_Z_Gat], axis = 0)
                    else:
                        delta_res = np.concatenate([res_X_Gat_C - solidInterfaceResidual0_X_Gat_C, res_Y_Gat_C - solidInterfaceResidual0_Y_Gat_C], axis=0)
                        delta_d = np.concatenate([solidInterfaceDisplacement_tilde_X_Gat - solidInterfaceDisplacement_tilde1_X_Gat, solidInterfaceDisplacement_tilde_Y_Gat - solidInterfaceDisplacement_tilde1_Y_Gat], axis = 0)
                    
                    self.Vk.insert(0, delta_res)
                    self.Wk.insert(0, delta_d)
                
                    V = np.vstack(self.Vk).T
                    W = np.vstack(self.Wk).T

                    V, W = self.filter(V, W)
                    U = np.transpose(W - np.dot(self.invJprev, V))
                    self.invJ = self.invJprev + np.linalg.lstsq(V.T, U, rcond=-1)[0].T

                if self.manager.nDim == 3:
                    Res = np.concatenate([res_X_Gat_C, res_Y_Gat_C, res_Z_Gat_C], axis=0)
                else:
                    Res = np.concatenate([res_X_Gat_C, res_Y_Gat_C], axis=0)

                if self.manager.nDim == 3:
                    delta_ds_loc = np.split(np.dot(self.invJ, -Res) + Res,3,axis=0)
                    
                    delta_ds_loc_X = delta_ds_loc[0]
                    delta_ds_loc_Y = delta_ds_loc[1]
                    delta_ds_loc_Z = delta_ds_loc[2]
                else:
                    delta_ds_loc = np.split(np.dot(self.invJ, -Res) + Res,2,axis=0)
                    
                    delta_ds_loc_X = delta_ds_loc[0]
                    delta_ds_loc_Y = delta_ds_loc[1]
                    delta_ds_loc_Z = np.zeros(ns)
                
                for iVertex in range(delta_ds_loc_X.shape[0]):
                    iGlobalVertex = self.manager.getGlobalIndex('solid', self.myid, iVertex)
                    delta_ds[iGlobalVertex] = [delta_ds_loc_X[iVertex], delta_ds_loc_Y[iVertex], delta_ds_loc_Z[iVertex]]
            
            # --- Go back to parallel run --- #
            mpiBarrier(self.mpiComm)
            delta_ds.assemble()
            self.interfaceInterpolator.solidInterfaceDisplacement += delta_ds
        
        if self.computeTangentMatrixBasedOnFirstIt:
            if self.FSIIter == 0:
                self.solidInterfaceResidual.copy(self.solidInterfaceResidual0)
                self.solidInterfaceDisplacement_tilde.copy(self.solidInterfaceDisplacement_tilde1)
        else:
            self.solidInterfaceResidual.copy(self.solidInterfaceResidual0)
            self.solidInterfaceDisplacement_tilde.copy(self.solidInterfaceDisplacement_tilde1)


    def relaxMultiVectorJacobian_CHT(self):

        if self.interfaceInterpolator.chtTransferMethod != 'FFTB': raise Exception('Only FFTB coupling is implemented with IQN-MVJ')

        d = self.interfaceInterpolator.getd()
        ns = self.interfaceInterpolator.getNs()
        delta_ds = FlexInterfaceData(ns+d, 1, self.mpiComm)

        # --- Initialize d_tilde for the construction of the WkCHT matrix -- #
        if self.myid in self.manager.getSolidInterfaceProcessors():
            localSolidInterfaceTemp = self.SolidSolver.getNodalTemperatures()
            for iVertex in range(self.manager.getNumberOfLocalSolidInterfaceNodes()):
                iGlobalVertex = self.manager.getGlobalIndex('solid', self.myid, iVertex)
                self.solidInterfaceTemperature_tilde[iGlobalVertex] = [localSolidInterfaceTemp[iVertex]]

        self.solidInterfaceTemperature_tilde.assemble()
        
        # --- Relax the solid position --- #
        if self.makeBGS_CHT:
            self.invJprevCHT = np.zeros((ns,ns))
            self.invJCHT = np.zeros((ns,ns))
            self.makeBGS_CHT = False
            self.relax_CHT()

        # --- Construct VkCHT and WkCHT matrices for the computation of the approximated tangent matrix --- #
        else:
            mpiPrint('\nCorrect solid interface temperature using IQN-MVJ method...\n', self.mpiComm)
            
            # --- Start gathering on root process --- #
            res_Gat = mpiGatherInterfaceData(self.solidTemperatureResidual, ns+d, self.mpiComm, 0)[0]
            solidInterfaceResidual0_Gat = mpiGatherInterfaceData(self.solidInterfaceResidual0_CHT, ns+d, self.mpiComm, 0)[0]
            solidInterfaceTemperature_tilde_Gat = mpiGatherInterfaceData(self.solidInterfaceTemperature_tilde, ns, self.mpiComm, 0)[0]
            solidInterfaceTemperature_tilde1_Gat = mpiGatherInterfaceData(self.solidInterfaceTemperature_tilde1, ns, self.mpiComm, 0)[0]

            # --- Copies for operating on, length=ns, not ns+d --- #
            if self.myid == 0:
                res_Gat_C = res_Gat[:ns]
                solidInterfaceResidual0_Gat_C = solidInterfaceResidual0_Gat[:ns]

                # --- Use J from the previous time step because VkCHT and WkCHT are empty --- #
                if self.FSIIter == 0:
                    self.invJCHT = self.invJprevCHT.copy()

                # -- VkCHT and WkCHT matrices are enriched only starting from the second iteration --- #
                else:
                    delta_res = res_Gat_C - solidInterfaceResidual0_Gat_C
                    delta_d = solidInterfaceTemperature_tilde_Gat - solidInterfaceTemperature_tilde1_Gat
                    
                    self.VkCHT.insert(0, delta_res)
                    self.WkCHT.insert(0, delta_d)
                
                    V = np.vstack(self.VkCHT).T
                    W = np.vstack(self.WkCHT).T

                    V, W = self.filter(V, W)
                    U = np.transpose(W - np.dot(self.invJprevCHT, V))
                    self.invJCHT = self.invJprevCHT + np.linalg.lstsq(V.T, U, rcond=-1)[0].T

                delta_ds_loc = np.dot(self.invJCHT, -res_Gat_C) + res_Gat_C

                for iVertex in range(delta_ds_loc.shape[0]):
                    iGlobalVertex = self.manager.getGlobalIndex('solid', self.myid, iVertex)
                    delta_ds[iGlobalVertex] = [delta_ds_loc[iVertex]]
            
            # --- Go back to parallel run --- #
            mpiBarrier(self.mpiComm)
            delta_ds.assemble()
            self.interfaceInterpolator.solidInterfaceTemperature += delta_ds
        
        if self.computeTangentMatrixBasedOnFirstIt:
            if self.FSIIter == 0:
                self.solidTemperatureResidual.copy(self.solidInterfaceResidual0_CHT)
                self.solidInterfaceTemperature_tilde.copy(self.solidInterfaceTemperature_tilde1)
        else:
            self.solidTemperatureResidual.copy(self.solidInterfaceResidual0_CHT)
            self.solidInterfaceTemperature_tilde.copy(self.solidInterfaceTemperature_tilde1)