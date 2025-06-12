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

RBMIntegrator.py
Python interface between the wrapper of NativeSolid and CUPyDO.
Authors D. THOMAS

'''

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import NativeSolid

import numpy as np
from ..utilities import titlePrint
from ..genericSolvers import SolidSolver, SolidAdjointSolver

# ----------------------------------------------------------------------
#  RBMI solver interface class
# ----------------------------------------------------------------------

class RBMI(SolidSolver):

    def __init__(self, p):


        self.NativeSolid = NativeSolid.NativeSolidSolver(p['csdFile'], True)
        self.regime = p['regime']

        self.interfaceID = self.NativeSolid.getFSIMarkerID()
        self.nNodes = self.NativeSolid.getNumberOfSolidInterfaceNodes(self.interfaceID)
        self.nHaloNode = 0
        self.nPhysicalNodes = self.NativeSolid.getNumberOfSolidInterfaceNodes(self.interfaceID)

        if self.regime == "harmonic":
            self.nInst = 2*self.NativeSolid.getNumberHarmonics()+1
        else:
            self.nInst = 1

        self.nDV = self.NativeSolid.getNumberDesignVariables()

        for iAlpha in range(int(self.nDV/2)): # Define location of centres of design variables
            x = 1-np.cos((iAlpha+1)*np.pi/(2*int(self.nDV/2)+1))
            self.NativeSolid.setDesignVariableCentre(x, 2*iAlpha)
            self.NativeSolid.setDesignVariableSide(1, 2*iAlpha)
            self.NativeSolid.setDesignVariableSide(-1, 2*iAlpha+1)
            self.NativeSolid.setDesignVariableCentre(x, 2*iAlpha+1)
        SolidSolver.__init__(self)

        self.__setInitialState()
        self.initRealTimeData()
        if not self.regime == 'harmonic':
            self.NativeSolid.updateSolution()

    def preprocessTimeIter(self, timeIter):

        self.NativeSolid.preprocessIteration(timeIter)

    def setInitialDisplacements(self):
        self.NativeSolid.setInitialDisplacements()
        self.__setInitialState()
        
    def run(self, t1, t2):


        if self.regime == 'unsteady':
            self.NativeSolid.timeIteration(t1, t2)
        elif self.regime == 'harmonic':
            self.NativeSolid.harmonicComputation()
        else:
            self.NativeSolid.staticComputation()

        self.__setCurrentState()
        return True

    def __setCurrentState(self):
        for jInst in range(self.nInst):
            if self.regime == "harmonic":
                self.NativeSolid.computeInterfacePosVel(False, jInst)
            for iVertex in range(self.nPhysicalNodes):
                self.nodalDisp_X[iVertex, jInst] = self.NativeSolid.getInterfaceNodeDispX(self.interfaceID, iVertex)
                self.nodalDisp_Y[iVertex, jInst] = self.NativeSolid.getInterfaceNodeDispY(self.interfaceID, iVertex)
                self.nodalDisp_Z[iVertex, jInst] = self.NativeSolid.getInterfaceNodeDispZ(self.interfaceID, iVertex)
                self.nodalVel_X[iVertex, jInst] = self.NativeSolid.getInterfaceNodeVelX(self.interfaceID, iVertex)
                self.nodalVel_Y[iVertex, jInst] = self.NativeSolid.getInterfaceNodeVelY(self.interfaceID, iVertex)
                self.nodalVel_Z[iVertex, jInst] = self.NativeSolid.getInterfaceNodeVelZ(self.interfaceID, iVertex)
    
    def __setInitialState(self):
        for jInst in range(self.nInst):
            self.NativeSolid.computeInterfacePosVel(False, jInst)
            for iVertex in range(self.nPhysicalNodes):
                self.nodalDisp_X[iVertex, jInst] = self.NativeSolid.getInterfaceNodeDispX(self.interfaceID, iVertex)
                self.nodalDisp_Y[iVertex, jInst] = self.NativeSolid.getInterfaceNodeDispY(self.interfaceID, iVertex)
                self.nodalDisp_Z[iVertex, jInst] = self.NativeSolid.getInterfaceNodeDispZ(self.interfaceID, iVertex)
                self.nodalVel_X[iVertex, jInst] = self.NativeSolid.getInterfaceNodeVelX(self.interfaceID, iVertex)
                self.nodalVel_Y[iVertex, jInst] = self.NativeSolid.getInterfaceNodeVelY(self.interfaceID, iVertex)
                self.nodalVel_Z[iVertex, jInst] = self.NativeSolid.getInterfaceNodeVelZ(self.interfaceID, iVertex)

    def getNodalInitialPositions(self):


        nodalInitialPos_X = np.zeros((self.nPhysicalNodes)) # initial position of the f/s interface
        nodalInitialPos_Y = np.zeros((self.nPhysicalNodes))
        nodalInitialPos_Z = np.zeros((self.nPhysicalNodes))

        for iVertex in range(self.nPhysicalNodes):
            nodalInitialPos_X[iVertex] = self.NativeSolid.getInterfaceNodePosX0(self.interfaceID, iVertex)
            nodalInitialPos_Y[iVertex] = self.NativeSolid.getInterfaceNodePosY0(self.interfaceID, iVertex)
            nodalInitialPos_Z[iVertex] = self.NativeSolid.getInterfaceNodePosZ0(self.interfaceID, iVertex)

        return (nodalInitialPos_X, nodalInitialPos_Y, nodalInitialPos_Z)

    def getNodalIndex(self, iVertex):


        self.NativeSolid.getInterfaceNodeGlobalIndex(self.interfaceID, iVertex)

    def applyNodalForce(self, load_X, load_Y, load_Z, dt, haloNodesLoads):


        for jInst in range(self.nInst):
            total_Y_load = 0.
            for iVertex in range(self.nPhysicalNodes):
                total_Y_load += load_Y[jInst][iVertex]
                self.NativeSolid.applyload(iVertex, jInst, load_X[jInst][iVertex], load_Y[jInst][iVertex], load_Z[jInst][iVertex])

            if self.regime == 'harmonic':
                self.NativeSolid.computeInterfacePosVel(False, jInst) # To update the position of the centre of rotation... inelegant
            self.NativeSolid.setGeneralisedForce(jInst)
            self.NativeSolid.setGeneralisedMoment(jInst)

    def getDeltaOmega(self):
        return self.NativeSolid.getDeltaOmega()
    
    def setOmegaHB(self, omega):
        if self.regime == 'harmonic':
            self.NativeSolid.setOmega(omega)

    def getNumberDesignVariables(self):
        return self.nDV

    def applyDesignVariables(self, alpha):
        for iAlpha in range(self.nDV):
            self.NativeSolid.setDesignVariableMagnitude(alpha[iAlpha], iAlpha)
        if alpha.size > self.nDV:
            self.NativeSolid.setDamping(0, alpha[-1])
        self.NativeSolid.applyDesignVariables()
        return
    
    def getObjectiveFunction(self):
        ObjFun = self.NativeSolid.getObjectiveFunction()
        return ObjFun

    def update(self):


        SolidSolver.update(self)

        self.NativeSolid.updateSolution()

    def saveRealTimeData(self, time, nFSIIter):


        self.NativeSolid.writeSolution(time, nFSIIter)

    def save(self):


        self.NativeSolid.saveSolution()

    def exit(self):


        print("***************************** Exit RBM Integrator *****************************")


class RBMIAdjoint(RBMI, SolidAdjointSolver):

    def __init__(self, p):
        self.NativeSolid = NativeSolid.NativeSolidSolver(p['csdFile'], True)
        self.regime = p['regime']

        self.interfaceID = self.NativeSolid.getFSIMarkerID()
        self.nNodes = self.NativeSolid.getNumberOfSolidInterfaceNodes(self.interfaceID)
        self.nHaloNode = 0
        self.nPhysicalNodes = self.NativeSolid.getNumberOfSolidInterfaceNodes(self.interfaceID)
        if self.regime == "harmonic":
            self.nInst = 2*self.NativeSolid.getNumberHarmonics()+1
        else:
            self.nInst = 1
        
        self.nDV = self.NativeSolid.getNumberDesignVariables()

        for iAlpha in range(int(self.nDV/2)): # Define location of centres of design variables
            x = 1-np.cos((iAlpha+1)*np.pi/(2*int(self.nDV/2)+1))
            self.NativeSolid.setDesignVariableCentre(x, 2*iAlpha)
            self.NativeSolid.setDesignVariableSide(1, 2*iAlpha)
            self.NativeSolid.setDesignVariableSide(-1, 2*iAlpha+1)
            self.NativeSolid.setDesignVariableCentre(x, 2*iAlpha+1)

        self.haloNodeList = {}

        # --- Create the array for external communication (displacement, velocity and velocity at the previous time step) --- #
        self.nodalDisp_X = np.zeros((self.nPhysicalNodes, self.nInst))
        self.nodalDisp_Y = np.zeros((self.nPhysicalNodes, self.nInst))
        self.nodalDisp_Z = np.zeros((self.nPhysicalNodes, self.nInst))
        self.nodalVel_X = np.zeros((self.nPhysicalNodes, self.nInst))
        self.nodalVel_Y = np.zeros((self.nPhysicalNodes, self.nInst))
        self.nodalVel_Z = np.zeros((self.nPhysicalNodes, self.nInst))

        self.nodalAdjDisp_X = np.zeros((self.nPhysicalNodes, self.nInst))
        self.nodalAdjDisp_Y = np.zeros((self.nPhysicalNodes, self.nInst))
        self.nodalAdjDisp_Z = np.zeros((self.nPhysicalNodes, self.nInst))

        self.nodalAdjLoad_X = np.zeros((self.nPhysicalNodes, self.nInst))
        self.nodalAdjLoad_Y = np.zeros((self.nPhysicalNodes, self.nInst))
        self.nodalAdjLoad_Z = np.zeros((self.nPhysicalNodes, self.nInst))

        RBMI._RBMI__setCurrentState(self)  # Using RBMI original __setCurrentState
        self.__setCurrentState()
        self.nodalVel_XNm1 = self.nodalVel_X.copy()
        self.nodalVel_YNm1 = self.nodalVel_Y.copy()
        self.nodalVel_ZNm1 = self.nodalVel_Z.copy()

        self.initRealTimeData()

    def applyNodalAdjointDisplacement(self, disp_adj_X, disp_adj_Y, disp_adj_Z, time, haloNodesDisplacements):
        for jInst in range(self.nInst):
            PhysicalIndex = 0
            for iVertex in range(self.nNodes):
                dispX = disp_adj_X[jInst][PhysicalIndex]
                dispY = disp_adj_Y[jInst][PhysicalIndex]
                dispZ = disp_adj_Z[jInst][PhysicalIndex]
                PhysicalIndex += 1
                self.NativeSolid.applyDisplacementAdjoint(iVertex, jInst, dispX, dispY, dispZ)
            self.NativeSolid.setTotalAdjointDisplacement(jInst)

    def applyFrequencyDerivative(self, omega_adj):
        self.NativeSolid.setFrequencyDerivative(omega_adj)

    def run(self, t1, t2):
        if self.regime == 'unsteady':
            self.NativeSolid.timeIteration(t1, t2)
        elif self.regime == 'harmonic':
            self.NativeSolid.harmonicComputation()
        else:
            self.NativeSolid.staticComputation()

        self.__setCurrentState()
        return True

    def __setCurrentState(self):
        for jInst in range(self.nInst):
            self.NativeSolid.computeInterfacePosVel(False, jInst)
            self.NativeSolid.computeInterfaceAdjointLoads(jInst)
            for iVertex in range(self.nPhysicalNodes):
                self.nodalAdjLoad_X[iVertex, jInst] = self.NativeSolid.getLoadDerivativeX(iVertex)
                self.nodalAdjLoad_Y[iVertex, jInst] = self.NativeSolid.getLoadDerivativeY(iVertex)
                self.nodalAdjLoad_Z[iVertex, jInst] = self.NativeSolid.getLoadDerivativeZ(iVertex)
    
    def getNodalAdjointLoads(self):
        return (self.nodalAdjLoad_X, self.nodalAdjLoad_Y, self.nodalAdjLoad_Z)

    def saveRealTimeData(self, time, nFSIIter):
        self.NativeSolid.writeAdjointSolution(time, nFSIIter)
    
    def getDesignVariableDerivative(self, iAlpha):
        return self.NativeSolid.getDesignVariableDerivative(iAlpha)
