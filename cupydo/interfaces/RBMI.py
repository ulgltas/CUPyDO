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

import math
import numpy as np
from ..utilities import titlePrint
from ..genericSolvers import SolidSolver

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

        SolidSolver.__init__(self)

        self.NativeSolid.setInitialDisplacements()
        self.__setCurrentState()
        self.initRealTimeData()

    def preprocessTimeIter(self, timeIter):

        self.NativeSolid.preprocessIteration(timeIter)
        
    def run(self, t1, t2):


        if self.regime == 'unsteady':
            self.NativeSolid.timeIteration(t1, t2)
        else:
            self.NativeSolid.staticComputation()

        self.__setCurrentState()
        return True

    def __setCurrentState(self):


        for iVertex in range(self.nPhysicalNodes):
            self.nodalDisp_X[iVertex] = self.NativeSolid.getInterfaceNodeDispX(self.interfaceID, iVertex)
            self.nodalDisp_Y[iVertex] = self.NativeSolid.getInterfaceNodeDispY(self.interfaceID, iVertex)
            self.nodalDisp_Z[iVertex] = self.NativeSolid.getInterfaceNodeDispZ(self.interfaceID, iVertex)
            self.nodalVel_X[iVertex] = self.NativeSolid.getInterfaceNodeVelX(self.interfaceID, iVertex)
            self.nodalVel_Y[iVertex] = self.NativeSolid.getInterfaceNodeVelY(self.interfaceID, iVertex)
            self.nodalVel_Z[iVertex] = self.NativeSolid.getInterfaceNodeVelZ(self.interfaceID, iVertex)

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


        for iVertex in range(self.nPhysicalNodes):
            self.NativeSolid.applyload(iVertex, load_X[iVertex], load_Y[iVertex], load_Z[iVertex])

        self.NativeSolid.setGeneralisedForce()
        self.NativeSolid.setGeneralisedMoment()

    def update(self):


        SolidSolver.update(self)

        self.NativeSolid.updateSolution()

    def saveRealTimeData(self, time, nFSIIter):


        self.NativeSolid.writeSolution(time, nFSIIter)

    def save(self):


        self.NativeSolid.saveSolution()

    def exit(self):


        titlePrint("Exit RBM Integrator")
