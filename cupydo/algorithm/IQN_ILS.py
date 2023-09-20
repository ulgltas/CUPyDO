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
import scipy as sp
import sys
import copy
import sys

from ..utilities import *
from ..interfaceData import FlexInterfaceData
from .BGS import AlgorithmBGSStaticRelax

np.set_printoptions(threshold=sys.maxsize)

# ----------------------------------------------------------------------
#    IQN ILS Algorithm class
# ----------------------------------------------------------------------

class AlgorithmIQN_ILS(AlgorithmBGSStaticRelax):
    def __init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, Criterion, p, mpiComm=None):

        AlgorithmBGSStaticRelax.__init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, Criterion, p, mpiComm)

        self.resetInternalVars()
        
        # --- Number of previous time steps used in the approximation of the tangent matrix --- #
        self.nbTimeToKeep = p['nSteps']

        # --- Option which allows to build the tangent matrix of a given time step using differences with respect to the first FSI iteration (delta_r_k = r_k+1 - r_0) instead of the previous iteration (delta_r_k = r_k+1 - r_k) --- #
        self.computeTangentMatrixBasedOnFirstIt = p['firstItTgtMat']

        # --- Tolerance and type of filtering : None, Degroote1, Degroote2, and Haelterman.
        self.tollQR = 1.0e-1
        self.qrFilter = 'Haelterman'
    
    def qrSolve(self, V, W, res):

        if self.qrFilter == None: # Classical least squares without QR filtering
            c = np.linalg.lstsq(V, -res, rcond=-1)[0]
        
        if self.qrFilter == 'Degroote1': # QR filtering as described by J. Degroote et al. Computers and Structures, 87, 793-801 (2009).
            Q, R = sp.linalg.qr(V, mode='economic')
            s = np.dot(np.transpose(Q), -res)
            toll = self.tollQR*sp.linalg.norm(R, 2)
            c = solve_upper_triangular_mod(R, s, toll)
        
        elif self.qrFilter == 'Degroote2': # QR filtering as described by J. Degroote et al. CMAME, 199, 2085-2098 (2010).
            Q, R, V, W = QRfiltering(V, W, self.tollQR)
            s = np.dot(np.transpose(Q), -res)
            c = np.linalg.solve(R, s)
        
        elif self.qrFilter == 'Haelterman': # 'Modified' QR filtering as described by R. Haelterman et al. Computers and Structures, 171, 9-17 (2016).
            Q, R, V, W = QRfiltering_mod(V, W, self.tollQR)
            s = np.dot(np.transpose(Q), -res)
            c = np.linalg.solve(R, s)
        
        else:
            raise NameError('IQN-ILS Algorithm: the QR filtering technique is unknown!')
        
        return c, W
    
    def resetInternalVars(self):

        # --- Indicate if a BGS iteration must be performed --- #
        if self.manager.mechanical:
            self.makeBGS = True
        if self.manager.thermal:
            self.makeBGS_CHT = True

            # --- Global V and W matrices for IQN-ILS algorithm, including information from previous time steps --- #
        if self.manager.mechanical:
            self.V = []
            self.W = []
        if self.manager.thermal:
            self.V_CHT = []
            self.W_CHT = []
        
    def fsiCoupling(self):
        """
        Interface Quasi Newton - Inverse Least Square (IQN-ILS) method for strong coupling FSI
        """

        mpiPrint("Enter IQN-ILS Strong Coupling FSI",self.mpiComm,titlePrint)

        self.FSIIter = 0
        self.FSIConv = False

        ns = self.interfaceInterpolator.getNs()
        d = self.interfaceInterpolator.getd()

        # --- Initialize all the quantities used in the IQN-ILS method --- #
        if self.manager.mechanical:
            self.solidInterfaceResidual0 = FlexInterfaceData(ns+d, 3, self.mpiComm)
            self.solidInterfaceDisplacement_tilde = FlexInterfaceData(ns, 3, self.mpiComm)
            self.solidInterfaceDisplacement_tilde1 = FlexInterfaceData(ns, 3, self.mpiComm)
            self.Vk_mat = np.zeros((self.manager.nDim*ns,1))
            self.Wk_mat = np.zeros((self.manager.nDim*ns,1))
            self.stack = 0
        if self.manager.thermal:
            self.solidInterfaceResidual0_CHT = FlexInterfaceData(ns+d, 1, self.mpiComm)
            self.solidInterfaceTemperature_tilde = FlexInterfaceData(ns, 1, self.mpiComm)
            self.solidInterfaceTemperature_tilde1 = FlexInterfaceData(ns, 1, self.mpiComm)
            self.Vk_matCHT = np.zeros((ns,1))
            self.Wk_matCHT = np.zeros((ns,1))
            self.stackCHT = 0

        if self.nbTimeToKeep > 0: # If information from previous time steps is re-used then Vk = V, Wk = W
            if self.manager.mechanical:
                self.Vk = copy.deepcopy(self.V)
                self.Wk = copy.deepcopy(self.W)
            if self.manager.thermal:
                self.VkCHT = copy.deepcopy(self.V_CHT)
                self.WkCHT = copy.deepcopy(self.W_CHT)
        else: # If information from previous time steps is not re-used then Vk and Wk are empty lists
            if self.manager.mechanical:
                self.Vk = []
                self.Wk = []
            if self.manager.thermal:
                self.VkCHT = []
                self.WkCHT = []

        while (self.FSIIter < self.nbFSIIterMax) and (self.FSIConv == False):
            mpiPrint("\n>>>> FSI iteration {} <<<<\n".format(self.FSIIter), self.mpiComm)

            if self.manager.mechanical:
                # --- Solid to fluid mechanical transfer --- #
                self.solidToFluidMechaTransfer()
                # --- Fluid mesh morphing --- #
                mpiPrint('\nPerforming mesh deformation...\n', self.mpiComm)
                self.meshDefTimer.start()
                self.FluidSolver.meshUpdate(self.step.timeIter)
                self.meshDefTimer.stop()
                self.meshDefTimer.cumul()

            if self.manager.thermal:
                # --- Solid to fluid thermal transfer --- #
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
            if not verif: return False

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
            if not verif: return False

            if self.manager.mechanical:
                # --- Compute and monitor the FSI residual --- #
                res = self.computeSolidInterfaceResidual()
                self.criterion.update(res)
                mpiPrint('\nFSI error value : {}\n'.format(self.criterion.epsilon), self.mpiComm)
                self.relaxInverseLeastSquare(res)

            if self.manager.thermal:
                # --- Compute and monitor the FSI residual --- #
                res = self.computeSolidInterfaceResidual_CHT()
                self.criterion.update_CHT(res)
                mpiPrint('\nCHT error value : {}\n'.format(self.criterion.epsilonCHT), self.mpiComm)
                self.relaxInverseLeastSquare_CHT(res)

            # --- Update the FSI iteration and history --- #
            if self.writeInFSIloop == True:
                self.writeRealTimeData()
            self.FSIIter += 1

            # --- Compute and monitor the FSI residual --- #
            if self.criterion.isVerified():
                mpiPrint("IQN-ILS is Converged",self.mpiComm,titlePrint)
                self.FSIConv = True
        
        # --- Empty the V and W containers if not converged --- #
        if not self.FSIConv: return False
        
        # ---  Add the current time data step to V and W, and remove out-of-range time steps --- #
        if self.manager.mechanical:
            if (self.nbTimeToKeep > 0) and (self.stack > 0):
                mpiPrint('\nUpdating V and W matrices...\n', self.mpiComm)
                self.V.insert(0, self.Vk_mat[:,0:self.stack].T)
                self.W.insert(0, self.Wk_mat[:,0:self.stack].T)
                while len(self.V) > self.nbTimeToKeep:
                    del self.V[-1]
                    del self.W[-1]

            # --- Start the next time step by a BGS since V and W are empty --- #
            if self.nbTimeToKeep == 0 or len(self.Vk) == 0:
                self.makeBGS = True

        # ---  Add the current time data step to V_CHT and W_CHT, and remove out-of-range time steps --- #
        if self.manager.thermal:
            if (self.nbTimeToKeep > 0) and (self.stackCHT > 0):
                mpiPrint('\nUpdating V_CHT and W_CHT matrices...\n', self.mpiComm)
                self.V_CHT.insert(0, self.Vk_matCHT[:,0:self.stackCHT].T)
                self.W_CHT.insert(0, self.Wk_matCHT[:,0:self.stackCHT].T)
                while len(self.V_CHT) > self.nbTimeToKeep:
                    del self.V_CHT[-1]
                    del self.W_CHT[-1]

            # --- Start the next time step by a BGS since V_CHT and W_CHT are empty --- #
            if self.nbTimeToKeep == 0 or len(self.VkCHT) == 0:
                self.makeBGS_CHT = True

        return True


    def relaxInverseLeastSquare(self, res):
            
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
            mpiPrint('\nProcessing interface displacement...\n', self.mpiComm)
            self.relaxSolidPosition()
            self.makeBGS = False

        # --- Construct Vk and Wk matrices for the computation of the approximated tangent matrix --- #
        else:
            mpiPrint('\nCorrect solid interface displacement using IQN-ILS method...\n', self.mpiComm)
            
            # --- Start gathering on root process --- #
            res_X_Gat, res_Y_Gat, res_Z_Gat = mpiGatherInterfaceData(res, ns+d, self.mpiComm, 0)
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

                # --- Vk and Wk matrices are enriched only starting from the second iteration of every FSI loop --- #
                if self.FSIIter > 0:
                    if self.manager.nDim == 3:
                        delta_res = np.concatenate([res_X_Gat_C - solidInterfaceResidual0_X_Gat_C, res_Y_Gat_C - solidInterfaceResidual0_Y_Gat_C, res_Z_Gat_C - solidInterfaceResidual0_Z_Gat_C], axis=0)
                        delta_d = np.concatenate([solidInterfaceDisplacement_tilde_X_Gat - solidInterfaceDisplacement_tilde1_X_Gat, solidInterfaceDisplacement_tilde_Y_Gat - solidInterfaceDisplacement_tilde1_Y_Gat, solidInterfaceDisplacement_tilde_Z_Gat - solidInterfaceDisplacement_tilde1_Z_Gat], axis = 0)
                    else:
                        delta_res = np.concatenate([res_X_Gat_C - solidInterfaceResidual0_X_Gat_C, res_Y_Gat_C - solidInterfaceResidual0_Y_Gat_C], axis=0)
                        delta_d = np.concatenate([solidInterfaceDisplacement_tilde_X_Gat - solidInterfaceDisplacement_tilde1_X_Gat, solidInterfaceDisplacement_tilde_Y_Gat - solidInterfaceDisplacement_tilde1_Y_Gat], axis = 0)
                    
                    self.Vk.insert(0, delta_res)
                    self.Wk.insert(0, delta_d)
                    self.stack += 1
                
                self.Vk_mat = np.vstack(self.Vk).T
                self.Wk_mat = np.vstack(self.Wk).T

                # --- Remove extra columns if number of iterations is larger than number of interface degrees of freedom --- #
                if (self.Vk_mat.shape[1] > self.manager.nDim*ns and self.qrFilter == 'Degroote1'):
                    mpiPrint('WARNING: IQN-ILS Algorithm using Degroote1 QR filter. Approximated stiffness matrix number of columns exceeds the number of degrees of freedom at FSI interface. Extra columns (the oldest ones!) are deleted for next iterations to avoid overdetermined problem!', self.mpiComm)
                    self.Vk_mat = np.delete(self.Vk_mat, np.s_[(self.manager.nDim*ns-self.Vk_mat.shape[1]):], 1)
                    self.Wk_mat = np.delete(self.Wk_mat, np.s_[(self.manager.nDim*ns-self.Wk_mat.shape[1]):], 1)
                
                dummy_V = self.Vk_mat.copy()
                dummy_W = self.Wk_mat.copy()
                
                if self.manager.nDim == 3:
                    dummy_Res = np.concatenate([res_X_Gat_C, res_Y_Gat_C, res_Z_Gat_C], axis=0)
                else:
                    dummy_Res = np.concatenate([res_X_Gat_C, res_Y_Gat_C], axis=0)
                
                c, dummy_W = self.qrSolve(dummy_V, dummy_W, dummy_Res)

                if self.manager.nDim == 3:
                    delta_ds_loc = np.split((np.dot(dummy_W,c).T + np.concatenate([res_X_Gat_C, res_Y_Gat_C, res_Z_Gat_C], axis=0)),3,axis=0)
                    
                    delta_ds_loc_X = delta_ds_loc[0]
                    delta_ds_loc_Y = delta_ds_loc[1]
                    delta_ds_loc_Z = delta_ds_loc[2]
                else:
                    delta_ds_loc = np.split((np.dot(dummy_W,c).T + np.concatenate([res_X_Gat_C, res_Y_Gat_C], axis=0)),2,axis=0)
                    
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
                res.copy(self.solidInterfaceResidual0)
                self.solidInterfaceDisplacement_tilde.copy(self.solidInterfaceDisplacement_tilde1)
        else:
            res.copy(self.solidInterfaceResidual0)
            self.solidInterfaceDisplacement_tilde.copy(self.solidInterfaceDisplacement_tilde1)


    def relaxInverseLeastSquare_CHT(self, res):

        if self.interfaceInterpolator.chtTransferMethod != 'FFTB': raise Exception('Only FFTB coupling is implemented with IQN-ILS')

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
            mpiPrint('\nProcessing interface temperature...\n', self.mpiComm)
            self.makeBGS_CHT = False
            self.relax_CHT()

        # --- Construct VkCHT and WkCHT matrices for the computation of the approximated tangent matrix --- #
        else:
            mpiPrint('\nCorrect solid interface temperature using IQN-ILS method...\n', self.mpiComm)
            
            # --- Start gathering on root process --- #
            res_Gat = mpiGatherInterfaceData(res, ns+d, self.mpiComm, 0)[0]
            solidInterfaceResidual0_Gat = mpiGatherInterfaceData(self.solidInterfaceResidual0_CHT, ns+d, self.mpiComm, 0)[0]
            solidInterfaceTemperature_tilde_Gat = mpiGatherInterfaceData(self.solidInterfaceTemperature_tilde, ns, self.mpiComm, 0)[0]
            solidInterfaceTemperature_tilde1_Gat = mpiGatherInterfaceData(self.solidInterfaceTemperature_tilde1, ns, self.mpiComm, 0)[0]

            # --- Copies for operating on, length=ns, not ns+d --- #
            if self.myid == 0:
                res_Gat_C = res_Gat[:ns]
                solidInterfaceResidual0_Gat_C = solidInterfaceResidual0_Gat[:ns]

                # --- VkCHT and WkCHT matrices are enriched only starting from the second iteration of every FSI loop --- #
                if self.FSIIter > 0:

                    delta_res = res_Gat_C - solidInterfaceResidual0_Gat_C
                    delta_d = solidInterfaceTemperature_tilde_Gat - solidInterfaceTemperature_tilde1_Gat
                    
                    self.VkCHT.insert(0, delta_res)
                    self.WkCHT.insert(0, delta_d)
                    self.stackCHT += 1
                
                self.Vk_matCHT = np.vstack(self.VkCHT).T
                self.Wk_matCHT = np.vstack(self.WkCHT).T
                
                # --- Remove extra columns if number of iterations is larger than number of interface degrees of freedom --- #
                if (self.Vk_matCHT.shape[1] > ns and self.qrFilter == 'Degroote1'):
                    mpiPrint('WARNING: IQN-ILS Algorithm using Degroote1 QR filter. Approximated stiffness matrix number of columns exceeds the number of degrees of freedom at FSI interface. Extra columns (the oldest ones!) are deleted for next iterations to avoid overdetermined problem!', self.mpiComm)
                    self.Vk_matCHT = np.delete(self.Vk_matCHT, np.s_[(ns-self.Vk_matCHT.shape[1]):], 1)
                    self.Wk_matCHT = np.delete(self.Wk_matCHT, np.s_[(ns-self.Wk_matCHT.shape[1]):], 1)
                
                dummy_V = self.Vk_matCHT.copy()
                dummy_W = self.Wk_matCHT.copy()
                c, dummy_W = self.qrSolve(dummy_V, dummy_W, res_Gat_C)
                delta_ds_loc = np.dot(dummy_W,c).T + res_Gat_C

                for iVertex in range(delta_ds_loc.shape[0]):
                    iGlobalVertex = self.manager.getGlobalIndex('solid', self.myid, iVertex)
                    delta_ds[iGlobalVertex] = [delta_ds_loc[iVertex]]
            
            # --- Go back to parallel run --- #
            mpiBarrier(self.mpiComm)
            delta_ds.assemble()
            self.interfaceInterpolator.solidInterfaceTemperature += delta_ds
        
        if self.computeTangentMatrixBasedOnFirstIt:
            if self.FSIIter == 0:
                res.copy(self.solidInterfaceResidual0_CHT)
                self.solidInterfaceTemperature_tilde.copy(self.solidInterfaceTemperature_tilde1)
        else:
            res.copy(self.solidInterfaceResidual0_CHT)
            self.solidInterfaceTemperature_tilde.copy(self.solidInterfaceTemperature_tilde1)