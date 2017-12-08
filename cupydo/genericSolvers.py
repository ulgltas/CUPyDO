#!/usr/bin/env python
# -*- coding: latin-1; -*-
#
# FSICoupler.py
# Main file (Python core) of CUPyDO.
# Authors : David THOMAS, Marco Lucio CERQUAGLIA, Romain BOMAN
#
# COPYRIGHT (C) University of Liège, 2017.

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

from math import *
import numpy as np
import scipy as sp
from scipy import spatial
import scipy.sparse.linalg as splinalg
import os, os.path, sys, string
import time as tm

import traceback

import socket, fnmatch
import fsi_pyutils

import copy

import ccupydo

np.set_printoptions(threshold=np.nan)

# global vars (underscore prevent them to be imported with "from module import *")
_theModule  = None
_theWDir    = None # workspace directory
_theWDirRoot = os.getcwd()  # base directory du workspace

# ----------------------------------------------------------------------
#  Generic solid solver class
# ----------------------------------------------------------------------

class SolidSolver:
    """
    Des.
    """

    def __init__(self):

        self.haloNodeList = {}

        # --- Create the array for external communication (displacement, velocity and velocity at the previous time step) --- #
        self.nodalDisp_X = np.zeros(self.nPhysicalNodes)
        self.nodalDisp_Y = np.zeros(self.nPhysicalNodes)
        self.nodalDisp_Z = np.zeros(self.nPhysicalNodes)
        self.nodalVel_X = np.zeros(self.nPhysicalNodes)
        self.nodalVel_Y = np.zeros(self.nPhysicalNodes)
        self.nodalVel_Z = np.zeros(self.nPhysicalNodes)
        self.nodalVel_XNm1 = np.zeros(self.nPhysicalNodes)
        self.nodalVel_YNm1 = np.zeros(self.nPhysicalNodes)
        self.nodalVel_ZNm1 = np.zeros(self.nPhysicalNodes)

        # --- Same for thermal coupling (heat fluxes and temperature) ---
        self.nodalHeatFlux_X = np.zeros(self.nPhysicalNodes)
        self.nodalHeatFlux_Y = np.zeros(self.nPhysicalNodes)
        self.nodalHeatFlux_Z = np.zeros(self.nPhysicalNodes)
        self.nodalTemperature = np.zeros(self.nPhysicalNodes)

    def setInitialDisplacements(self):
        return

    def preprocessTimeIter(self, timeIter):
        return

    def run(self):
        return

    def __setCurrentState(self):
        return

    def getNodalDisplacements(self):
        """
        Des.
        """

        return (self.nodalDisp_X, self.nodalDisp_Y, self.nodalDisp_Z)

    def getNodalHeatFluxes(self):
        """
        Des.
        """

        return (self.nodalHeatFlux_X, self.nodalHeatFlux_Y, self.nodalHeatFlux_Z)

    def getNodalTemperatures(self):
        """
        Des.
        """

        return self.nodalTemperature

    def getNodalInitialPositions(self):
        return

    def getNodalVelocity(self):
        """
        des.
        """

        return (self.nodalVel_X, self.nodalVel_Y, self.nodalVel_Z)

    def getNodalVelocityNm1(self):
        """
        Des.
        """

        return (self.nodalVel_XNm1, self.nodalVel_YNm1, self.nodalVel_ZNm1)

    def getNodalIndex(self, iVertex):
        return

    def fakeFluidSolver(self, time):
        return

    def applyNodalLoads(self, load_X, load_Y, load_Z, time):
        return

    def applyNodalTemperatures(self, Temperature, time):
        return

    def applyNodalNormalHeatFluxes(self, NormalHeatFlux, val_time):
        return

    def applyNodalHeatFluxes(self, HeatFlux_X, HeatFlux_Y, HeatFlux_Z, time):
        return

    def update(self):

        self.nodalVel_XNm1 = self.nodalVel_X.copy()
        self.nodalVel_YNm1 = self.nodalVel_Y.copy()
        self.nodalVel_ZNm1 = self.nodalVel_Z.copy()

    def bgsUpdate(self):
        return

    def save(self):
        return

    def initRealTimeData(self):
        return

    def saveRealTimeData(self, time, nFSIIter):
        return

    def printRealTimeData(self, time, nFSIIter):
        return

    def remeshing(self):
        return

    def exit(self):
        return

# ----------------------------------------------------------------------
#  Generic fluid solver class
# ----------------------------------------------------------------------

class FluidSolver:
    """
    Des.
    """

    def __init__(self):

        self.haloNodeList = {}

        self.nodalLoad_X = np.zeros((self.nPhysicalNodes))
        self.nodalLoad_Y = np.zeros((self.nPhysicalNodes))
        self.nodalLoad_Z = np.zeros((self.nPhysicalNodes))

        self.nodalTemperature = np.zeros((self.nPhysicalNodes))

        self.nodalNormalHeatFlux = np.zeros(self.nPhysicalNodes)

        self.nodalHeatFlux_X = np.zeros((self.nPhysicalNodes))
        self.nodalHeatFlux_Y = np.zeros((self.nPhysicalNodes))
        self.nodalHeatFlux_Z = np.zeros((self.nPhysicalNodes))

        self.QWallInit = 0
        self.TWallInit = 288.0

    def setInitialMeshDeformation(self):
        return

    def setInitialInterfaceHeatFlux(self):
        return

    def setInitialInterfaceTemperature(self):
        return

    def preprocessTimeIter(self, timeIter):
        return

    def run(self):
        return

    def getNodalIndex(self, iVertex):
        return

    def fakeSolidSolver(self, time):
        return

    def getNodalLoads(self):

        return (self.nodalLoad_X, self.nodalLoad_Y, self.nodalLoad_Z)

    def getNodalInitialPositions(self):
        return

    def getNodalTemperatures(self):
        return self.nodalTemperature

    def getNodalNormalHeatFlux(self):
        return self.nodalNormalHeatFlux

    def getNodalHeatFluxes(self):
        return (self.nodalHeatFlux_X, self.nodalHeatFlux_Y, self.nodalHeatFlux_Z)

    def applyNodalDisplacements(self, dx, dy, dz, dx_nM1, dy_nM1, dz_nM1, haloNodesDisplacements,time):
        return

    def applyNodalHeatFluxes(self, HF_X, HF_Y, HF_Z, time):
        return

    def applyNodalTemperatures(self, Temperature, time):
        return

    def update(self, dt):
        return

    def bgsUpdate(self):
        return

    def save(self, nt):
        return

    def initRealTimeData(self):
        return

    def saveRealTimeData(self, time, nFSIIter):
        return

    def printRealTimeData(self, time, nFSIIter):
        return

    def remeshing(self):
        return

    def meshUpdate(self, nt):
        return

    def exit(self):
        return
