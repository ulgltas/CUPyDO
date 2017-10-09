#!/usr/bin/env python
# -*- coding: latin-1; -*-

# RBMIntegratorInterface.py
# Python interface between the wrapper of NativeSolid and CUPyDO.
# Authors D. THOMAS
#
# COPYRIGHT (C) University of Liège, 2017.

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import NativeSolid

import math
import numpy as np
from FSICoupler import SolidSolver

# ----------------------------------------------------------------------
#  RBMIntegrator class
# ----------------------------------------------------------------------

class RBMIntegrator(SolidSolver):

    def __init__(self, confFile, computationType):
        """
        Des.
        """

        self.NativeSolid = NativeSolid.NativeSolidSolver(confFile, True)
        self.computationType = computationType

        self.interfaceID = self.NativeSolid.getFSIMarkerID()
        self.nNodes = self.NativeSolid.getNumberOfSolidInterfaceNodes(self.interfaceID)
        self.nHaloNode = 0
        self.nPhysicalNodes = self.NativeSolid.getNumberOfSolidInterfaceNodes(self.interfaceID)

        SolidSolver.__init__(self)

        self.__setCurrentState()
        self.nodalVel_XNm1 = self.nodalVel_X.copy()
        self.nodalVel_YNm1 = self.nodalVel_Y.copy()
        self.nodalVel_ZNm1 = self.nodalVel_Z.copy()

        self.initRealTimeData()

    def preprocessTimeIter(self, timeIter):
        """
        Des.
        """

        self.NativeSolid.preprocessIteration(timeIter)

    def setInitialDisplacements(self):
        """
        Des.
        """

        self.NativeSolid.setInitialDisplacements()
        self.__setCurrentState()
        
    def run(self, t1, t2):
        """
        Des.
        """

        if self.computationType == 'unsteady':
            self.NativeSolid.timeIteration(t1, t2)
        else:
            self.NativeSolid.staticComputation()

        self.__setCurrentState()

    def __setCurrentState(self):
        """
        Des.
        """

        for iVertex in range(self.nPhysicalNodes):
            self.nodalDisp_X[iVertex] = self.NativeSolid.getInterfaceNodeDispX(self.interfaceID, iVertex)
            self.nodalDisp_Y[iVertex] = self.NativeSolid.getInterfaceNodeDispY(self.interfaceID, iVertex)
            self.nodalDisp_Z[iVertex] = self.NativeSolid.getInterfaceNodeDispZ(self.interfaceID, iVertex)
            self.nodalVel_X[iVertex] = self.NativeSolid.getInterfaceNodeVelX(self.interfaceID, iVertex)
            self.nodalVel_Y[iVertex] = self.NativeSolid.getInterfaceNodeVelY(self.interfaceID, iVertex)
            self.nodalVel_Z[iVertex] = self.NativeSolid.getInterfaceNodeVelZ(self.interfaceID, iVertex)

    def getNodalInitialPositions(self):
        """
        des.
        """

        nodalInitialPos_X = np.zeros((self.nPhysicalNodes)) # initial position of the f/s interface
        nodalInitialPos_Y = np.zeros((self.nPhysicalNodes))
        nodalInitialPos_Z = np.zeros((self.nPhysicalNodes))

        for iVertex in range(self.nPhysicalNodes):
            nodalInitialPos_X[iVertex] = self.NativeSolid.getInterfaceNodePosX0(self.interfaceID, iVertex)
            nodalInitialPos_Y[iVertex] = self.NativeSolid.getInterfaceNodePosY0(self.interfaceID, iVertex)
            nodalInitialPos_Z[iVertex] = self.NativeSolid.getInterfaceNodePosZ0(self.interfaceID, iVertex)

        return (nodalInitialPos_X, nodalInitialPos_Y, nodalInitialPos_Z)

    def getNodalIndex(self, iVertex):
        """
        Des.
        """

        self.NativeSolid.getInterfaceNodeGlobalIndex(self.interfaceID, iVertex)

    def applyNodalLoads(self, load_X, load_Y, load_Z, time):
        """
        Des.
        """

        for iVertex in range(self.nPhysicalNodes):
            self.NativeSolid.applyload(iVertex, load_X[iVertex], load_Y[iVertex], load_Z[iVertex])

        self.NativeSolid.setGeneralisedForce()
        self.NativeSolid.setGeneralisedMoment()

    def update(self):
        """
        Des.
        """

        SolidSolver.update(self)

        self.NativeSolid.updateSolution()

    def initRealTimeData(self):
        """
        Des.
        """

        solFile = open('NodalDisplacement.ascii', "w")
        solFile.write("Time\tnIter\tValue\n")
        solFile.close()

    def saveRealTimeData(self, time, nFSIIter):
        """
        Des.
        """

        solFile = open('NodalDisplacement.ascii', "a")
        solFile.write(str(time) + '\t' + str(nFSIIter) + '\t' + str(0.0) + '\n')

    def save(self):
        """
        Des.
        """

        self.NativeSolid.saveSolution()

    def exit(self):
        """
        Des.
        """

        print("***************************** Exit RBM Integrator *****************************")
