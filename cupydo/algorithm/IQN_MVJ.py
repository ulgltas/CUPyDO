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
    """
    Des.
    """

    def __init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, Criterion, nbFSIIterMax, deltaT, totTime, dtSave, omegaBoundList= [1.0, 1.0], computeTangentMatrixBasedOnFirstIt = False, mpiComm=None):
        """
        Des.
        """

        AlgorithmBGSStaticRelax.__init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, Criterion, nbFSIIterMax, deltaT, totTime, dtSave, omegaBoundList, mpiComm)

        # --- Indicate if a BGS iteration must be performed --- #
        self.makeBGS = True
        self.hasInvJ = False

        # --- Internal variables for convcergence check --- #
        self.maxNbOfItReached = False
        self.convergenceReachedInOneIt = False

        # --- Option which allows to build the tangent matrix of a given time step using differences with respect to the first FSI iteration (delta_r_k = r_k+1 - r_0) instead of the previous iteration (delta_r_k = r_k+1 - r_k) --- #
        self.computeTangentMatrixBasedOnFirstIt = computeTangentMatrixBasedOnFirstIt
        
        # --- Global V and W matrices for IQN-MVJ algorithm --- #
        self.V = []
        self.W = []

        # --- Curent inverse approximate Jacobian (J^-1 - I) --- #
        ns = self.interfaceInterpolator.getNs()
    
    def fsiCoupling(self):
        """
        Interface Quasi Newton - Multi Vector Jacobian (IQN-MVJ) method for strong coupling FSI
        """

        mpiPrint("Enter IQN-MNV Strong Coupling FSI",self.mpiComm,titlePrint)

        self.FSIIter = 0
        self.FSIConv = False
        self.errValue = 1.0
        self.errValue_CHT = 1e6 # Just for compatibility
        
        ns = self.interfaceInterpolator.getNs()
        d = self.interfaceInterpolator.getd()

        # --- Initialize all the quantities used in the IQN-MVJ method --- #
        res = FlexInterfaceData(ns, 3, self.mpiComm)
        solidInterfaceResidual0 = FlexInterfaceData(ns+d, 3, self.mpiComm)

        solidInterfaceDisplacement_tilde = FlexInterfaceData(ns, 3, self.mpiComm)
        solidInterfaceDisplacement_tilde1 = FlexInterfaceData(ns, 3, self.mpiComm)

        delta_ds = FlexInterfaceData(ns+d, 3, self.mpiComm)

        delta_ds_loc_X = np.zeros(0)
        delta_ds_loc_Y = np.zeros(0)
        delta_ds_loc_Z = np.zeros(0)

        # --- Initialises the previous approximate inverse Jacobian and V, W --- #
        Vk = []
        Wk = []
        nIt = 0

        while (self.FSIIter < self.nbFSIIterMax):
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

            # --- The fluid solver failed if verif is false --- #
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

            # --- The solid solver failed if verif is false --- #
            try: solidProc = int(self.manager.getSolidSolverProcessors())
            except: raise Exception('Only one solid solver process is supported yet')
            verif = mpiScatter(verif, self.mpiComm, solidProc)
            self.solidHasRun = True
            if not verif: return False

            # --- Compute and monitor the FSI residual --- #
            res = self.computeSolidInterfaceResidual()
            self.errValue = self.criterion.update(res)
            mpiPrint('\nFSI error value : {}\n'.format(self.errValue), self.mpiComm)
            self.FSIConv = self.criterion.isVerified(self.errValue)

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
                self.invJprev = np.zeros((self.manager.nDim*ns,self.manager.nDim*ns))
                self.relaxSolidPosition()
                self.makeBGS = False

            else:
                # --- Construct Vk and Wk matrices for the computation of the approximated tangent matrix --- #
                mpiPrint('\nCorrect solid interface displacements using IQN-MVJ method...\n', self.mpiComm)
                
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

                    if self.FSIIter == 0: # Use J^-1 from the previous time step because Vk and Wk are empty
                        
                        self.invJ = self.invJprev.copy()
                        self.hasInvJ = True

                    else: # Vk and Wk matrices are enriched only starting from the second iteration of every FSI loop
                        if self.manager.nDim == 3:
                            delta_res = np.concatenate([res_X_Gat_C - solidInterfaceResidual0_X_Gat_C, res_Y_Gat_C - solidInterfaceResidual0_Y_Gat_C, res_Z_Gat_C - solidInterfaceResidual0_Z_Gat_C], axis=0)
                            delta_d = np.concatenate([solidInterfaceDisplacement_tilde_X_Gat - solidInterfaceDisplacement_tilde1_X_Gat, solidInterfaceDisplacement_tilde_Y_Gat - solidInterfaceDisplacement_tilde1_Y_Gat, solidInterfaceDisplacement_tilde_Z_Gat - solidInterfaceDisplacement_tilde1_Z_Gat], axis = 0)
                        else:
                            delta_res = np.concatenate([res_X_Gat_C - solidInterfaceResidual0_X_Gat_C, res_Y_Gat_C - solidInterfaceResidual0_Y_Gat_C], axis=0)
                            delta_d = np.concatenate([solidInterfaceDisplacement_tilde_X_Gat - solidInterfaceDisplacement_tilde1_X_Gat, solidInterfaceDisplacement_tilde_Y_Gat - solidInterfaceDisplacement_tilde1_Y_Gat], axis = 0)
                        
                        Vk.insert(0, delta_res)
                        Wk.insert(0, delta_d)
                        
                        nIt+=1
                    
                        Vk_mat = np.vstack(Vk).T
                        Wk_mat = np.vstack(Wk).T
                        
                        X = np.transpose(Wk_mat-np.dot(self.invJprev,Vk_mat))
                        self.invJ = self.invJprev+np.linalg.lstsq(Vk_mat.T,X,rcond=-1)[0].T
                        self.hasInvJ = True

                    if self.manager.nDim == 3:
                        Res = np.concatenate([res_X_Gat_C, res_Y_Gat_C, res_Z_Gat_C], axis=0)
                    else:
                        Res = np.concatenate([res_X_Gat_C, res_Y_Gat_C], axis=0)

                    if self.manager.nDim == 3:
                        delta_ds_loc = np.split(np.dot(self.invJ,-Res)+Res,3,axis=0)
                        
                        delta_ds_loc_X = delta_ds_loc[0]
                        delta_ds_loc_Y = delta_ds_loc[1]
                        delta_ds_loc_Z = delta_ds_loc[2]
                    else:
                        delta_ds_loc = np.split(np.dot(self.invJ,-Res)+Res,2,axis=0)
                        
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
            self.FSIIter += 1
            if self.criterion.isVerified(self.errValue, self.errValue_CHT):

                mpiPrint("MVJ is Converged",self.mpiComm,titlePrint)
                if self.hasInvJ: self.invJprev = np.copy(self.invJ)
                return True
        
        self.makeBGS = True
        return False
