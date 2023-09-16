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

        # --- Number of previous time steps used in the approximation of the tangent matrix --- #
        self.nbTimeToKeep = p['nSteps']

        # --- Option which allows to build the tangent matrix of a given time step using differences with respect to the first FSI iteration (delta_r_k = r_k+1 - r_0) instead of the previous iteration (delta_r_k = r_k+1 - r_k) --- #
        self.computeTangentMatrixBasedOnFirstIt = p['firstItTgtMat']

        # --- Tolerance and type of filtering : None, Degroote1, Degroote2, and Haelterman.
        self.tollQR = 1.0e-1
        self.qrFilter = 'Haelterman'

        # --- Indicate if a BGS iteration must be performed --- #
        self.makeBGS = True
        
        # --- Global V and W matrices for IQN-ILS algorithm, including information from previous time steps --- #
        self.V = []
        self.W = []
    
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

        self.makeBGS = True
        self.V = []
        self.W = []
    
    def fsiCoupling(self):
        """
        Interface Quasi Newton - Inverse Least Square (IQN-ILS) method for strong coupling FSI
        """

        mpiPrint("Enter IQN-ILS Strong Coupling FSI",self.mpiComm,titlePrint)

        stack = 0
        self.FSIIter = 0
        self.FSIConv = False
        self.errValue = 1.0
        self.errValue_CHT = 1e6 # Just for compatibility. CHT not yet implemented for the IQN-ILS algorithm.
        
        ns = self.interfaceInterpolator.getNs()
        d = self.interfaceInterpolator.getd()

        # --- Initialize all the quantities used in the IQN-ILS method --- #
        res = FlexInterfaceData(ns, 3, self.mpiComm)
        solidInterfaceResidual0 = FlexInterfaceData(ns+d, 3, self.mpiComm)

        solidInterfaceDisplacement_tilde = FlexInterfaceData(ns, 3, self.mpiComm)
        solidInterfaceDisplacement_tilde1 = FlexInterfaceData(ns, 3, self.mpiComm)

        delta_ds = FlexInterfaceData(ns+d, 3, self.mpiComm)

        Vk_mat = np.zeros((self.manager.nDim*ns,1))
        Wk_mat = np.zeros((self.manager.nDim*ns,1))

        delta_ds_loc_X = np.zeros(0)
        delta_ds_loc_Y = np.zeros(0)
        delta_ds_loc_Z = np.zeros(0)

        if self.nbTimeToKeep > 0: # If information from previous time steps is re-used then Vk = V, Wk = W
            Vk = copy.deepcopy(self.V)
            Wk = copy.deepcopy(self.W)
        else: # If information from previous time steps is not re-used then Vk and Wk are empty lists
            Vk = []
            Wk = []

        while (self.FSIIter < self.nbFSIIterMax) and (self.FSIConv == False):
            mpiPrint("\n>>>> FSI iteration {} <<<<\n".format(self.FSIIter), self.mpiComm)

            # --- Solid to fluid mechanical transfer --- #
            self.solidToFluidMechaTransfer()

            # --- Fluid mesh morphing --- #
            mpiPrint('\nPerforming mesh deformation...\n', self.mpiComm)
            self.meshDefTimer.start()
            self.FluidSolver.meshUpdate(self.step.timeIter)
            self.meshDefTimer.stop()
            self.meshDefTimer.cumul()

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
            mpiPrint('\nProcessing interface fluid loads...\n', self.mpiComm)
            self.fluidToSolidMechaTransfer()
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

            # --- Compute and monitor the FSI residual --- #
            res = self.computeSolidInterfaceResidual()
            self.errValue = self.criterion.update(res)
            mpiPrint('\nFSI error value : {}\n'.format(self.errValue), self.mpiComm)

            # --- Initialize d_tilde for the construction of the Wk matrix -- #
            if self.myid in self.manager.getSolidInterfaceProcessors():
                localSolidInterfaceDisp_X, localSolidInterfaceDisp_Y, localSolidInterfaceDisp_Z = self.SolidSolver.getNodalDisplacements()
                for iVertex in range(self.manager.getNumberOfLocalSolidInterfaceNodes()):
                    iGlobalVertex = self.manager.getGlobalIndex('solid', self.myid, iVertex)
                    solidInterfaceDisplacement_tilde[iGlobalVertex] = [localSolidInterfaceDisp_X[iVertex], localSolidInterfaceDisp_Y[iVertex], localSolidInterfaceDisp_Z[iVertex]]

            solidInterfaceDisplacement_tilde.assemble()
            
            if self.makeBGS:
                # --- Relax the solid position --- #
                mpiPrint('\nProcessing interface displacements...\n', self.mpiComm)
                self.relaxSolidPosition()
                self.makeBGS = False

            else:
                
                # --- Construct Vk and Wk matrices for the computation of the approximated tangent matrix --- #
                mpiPrint('\nCorrect solid interface displacements using IQN-ILS method...\n', self.mpiComm)
                
                # --- Start gathering on root process --- #
                res_X_Gat, res_Y_Gat, res_Z_Gat = mpiGatherInterfaceData(res, ns+d, self.mpiComm, 0)
                solidInterfaceResidual0_X_Gat, solidInterfaceResidual0_Y_Gat, solidInterfaceResidual0_Z_Gat = mpiGatherInterfaceData(solidInterfaceResidual0, ns+d, self.mpiComm, 0)
                solidInterfaceDisplacement_tilde_X_Gat, solidInterfaceDisplacement_tilde_Y_Gat, solidInterfaceDisplacement_tilde_Z_Gat = mpiGatherInterfaceData(solidInterfaceDisplacement_tilde, ns, self.mpiComm, 0)
                solidInterfaceDisplacement_tilde1_X_Gat, solidInterfaceDisplacement_tilde1_Y_Gat, solidInterfaceDisplacement_tilde1_Z_Gat = mpiGatherInterfaceData(solidInterfaceDisplacement_tilde1, ns, self.mpiComm, 0)
                if self.myid == 0:
                    res_X_Gat_C = res_X_Gat[:ns] # Copies for operating on, length=ns, not ns+d
                    res_Y_Gat_C = res_Y_Gat[:ns]
                    res_Z_Gat_C = res_Z_Gat[:ns]
                    solidInterfaceResidual0_X_Gat_C = solidInterfaceResidual0_X_Gat[:ns] # Copies for operating on, length=ns, not ns+d
                    solidInterfaceResidual0_Y_Gat_C = solidInterfaceResidual0_Y_Gat[:ns]
                    solidInterfaceResidual0_Z_Gat_C = solidInterfaceResidual0_Z_Gat[:ns]
                    if self.FSIIter > 0: # Either information from previous time steps is re-used or not, Vk and Wk matrices are enriched only starting from the second iteration of every FSI loop
                        if self.manager.nDim == 3:
                            delta_res = np.concatenate([res_X_Gat_C - solidInterfaceResidual0_X_Gat_C, res_Y_Gat_C - solidInterfaceResidual0_Y_Gat_C, res_Z_Gat_C - solidInterfaceResidual0_Z_Gat_C], axis=0)
                            delta_d = np.concatenate([solidInterfaceDisplacement_tilde_X_Gat - solidInterfaceDisplacement_tilde1_X_Gat, solidInterfaceDisplacement_tilde_Y_Gat - solidInterfaceDisplacement_tilde1_Y_Gat, solidInterfaceDisplacement_tilde_Z_Gat - solidInterfaceDisplacement_tilde1_Z_Gat], axis = 0)
                        else:
                            delta_res = np.concatenate([res_X_Gat_C - solidInterfaceResidual0_X_Gat_C, res_Y_Gat_C - solidInterfaceResidual0_Y_Gat_C], axis=0)
                            delta_d = np.concatenate([solidInterfaceDisplacement_tilde_X_Gat - solidInterfaceDisplacement_tilde1_X_Gat, solidInterfaceDisplacement_tilde_Y_Gat - solidInterfaceDisplacement_tilde1_Y_Gat], axis = 0)
                        
                        Vk.insert(0, delta_res)
                        Wk.insert(0, delta_d)
                        stack += 1
                    
                    Vk_mat = np.vstack(Vk).T
                    Wk_mat = np.vstack(Wk).T
                    
                    if (Vk_mat.shape[1] > self.manager.nDim*ns and self.qrFilter == 'Degroote1'): # Remove extra columns if number of iterations (i.e. columns of Vk and Wk) is larger than number of interface degrees of freedom 
                        mpiPrint('WARNING: IQN-ILS Algorithm using \'Degroote1\' QR filter. Approximated stiffness matrix number of columns exceeds the number of degrees of freedom at FSI interface. Extra columns (the oldest ones!) are deleted for next iterations to avoid overdetermined problem!', self.mpiComm)
                        Vk_mat = np.delete(Vk_mat, np.s_[(self.manager.nDim*ns-Vk_mat.shape[1]):], 1)
                        Wk_mat = np.delete(Wk_mat, np.s_[(self.manager.nDim*ns-Wk_mat.shape[1]):], 1)
                    
                    dummy_V = Vk_mat.copy()
                    dummy_W = Wk_mat.copy()
                    
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
                    res.copy(solidInterfaceResidual0)
                    solidInterfaceDisplacement_tilde.copy(solidInterfaceDisplacement_tilde1)
            else:
                res.copy(solidInterfaceResidual0)
                solidInterfaceDisplacement_tilde.copy(solidInterfaceDisplacement_tilde1)

            # --- Update the FSI iteration and history --- #
            if self.writeInFSIloop == True:
                self.writeRealTimeData()
            self.FSIIter += 1

            # --- Compute and monitor the FSI residual --- #
            if self.criterion.isVerified(self.errValue):
                mpiPrint("IQN-ILS is Converged",self.mpiComm,titlePrint)
                self.FSIConv = True
        
        # --- Empty the V and W containers if not converged --- #
        if not self.FSIConv: return False
        
        # Add the current time data step to V and W, and remove out-of-range time steps
        if (self.nbTimeToKeep > 0) and (stack > 0):
            
            mpiPrint('\nUpdating V and W matrices...\n', self.mpiComm)
            self.V.insert(0, Vk_mat[:,0:stack].T)
            self.W.insert(0, Wk_mat[:,0:stack].T)
            
            while len(self.V) > self.nbTimeToKeep:
                del self.V[-1]
                del self.W[-1]

        if self.nbTimeToKeep == 0:
            # --- Start the next time step by a BGS since V and W are empty --- #
            self.makeBGS = True
        return True
    