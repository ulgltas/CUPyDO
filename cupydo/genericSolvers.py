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

genericSolvers.py
Generic solver classes.
Authors : David THOMAS, Marco Lucio CERQUAGLIA, Romain BOMAN

'''

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import numpy as np
import sys

np.set_printoptions(threshold=sys.maxsize)

# ----------------------------------------------------------------------
#  Generic solid solver class
# ----------------------------------------------------------------------

class SolidSolver(object):
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
        return (self.nodalDisp_X, self.nodalDisp_Y, self.nodalDisp_Z)

    def getNodalHeatFluxes(self):
        return (self.nodalHeatFlux_X, self.nodalHeatFlux_Y, self.nodalHeatFlux_Z)

    def getNodalTemperatures(self):
        return self.nodalTemperature

    def getNodalInitialPositions(self):
        return

    def getNodalVelocity(self):
        return (self.nodalVel_X, self.nodalVel_Y, self.nodalVel_Z)

    def getNodalVelocityNm1(self):
        return (self.nodalVel_XNm1, self.nodalVel_YNm1, self.nodalVel_ZNm1)

    def getNodalIndex(self, iVertex):
        return

    def fakeFluidSolver(self, dt):
        return

    def applyNodalLoads(self, load_X, load_Y, load_Z, dt, haloNodesLoads = {}):
        return

    def applyNodalTemperatures(self, Temperature, dt):
        return

    def applyNodalNormalHeatFluxes(self, NormalHeatFlux, dt):
        return

    def applyNodalHeatFluxes(self, HeatFlux_X, HeatFlux_Y, HeatFlux_Z, dt):
        return

    def update(self):

        self.nodalVel_XNm1 = self.nodalVel_X.copy()
        self.nodalVel_YNm1 = self.nodalVel_Y.copy()
        self.nodalVel_ZNm1 = self.nodalVel_Z.copy()

    def steadyUpdate(self):
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

class FluidSolver(object):
    def __init__(self, p):

        self.haloNodeList = {}

        if p['interpType'] == 'conservative':
            self.nodalLoad_X = np.zeros((self.nPhysicalNodes))
            self.nodalLoad_Y = np.zeros((self.nPhysicalNodes))
            self.nodalLoad_Z = np.zeros((self.nPhysicalNodes))

        else:
            self.nodalLoad_XX = np.zeros((self.nPhysicalNodes))
            self.nodalLoad_YY = np.zeros((self.nPhysicalNodes))
            self.nodalLoad_ZZ = np.zeros((self.nPhysicalNodes))
            self.nodalLoad_XY = np.zeros((self.nPhysicalNodes))
            self.nodalLoad_XZ = np.zeros((self.nPhysicalNodes))
            self.nodalLoad_YZ = np.zeros((self.nPhysicalNodes))

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

    def run(self, t1, t2):
        return

    def getNodalIndex(self, iVertex):
        return

    def fakeSolidSolver(self, dt):
        return

    def getNodalLoads(self):

        return (self.nodalLoad_X, self.nodalLoad_Y, self.nodalLoad_Z)

    def getNodalStress(self):

        return (self.nodalLoad_XX, self.nodalLoad_YY, self.nodalLoad_ZZ, self.nodalLoad_XY, self.nodalLoad_XZ, self.nodalLoad_YZ)
    
    def getNodalInitialPositions(self):
        return

    def getNodalTemperatures(self):
        return self.nodalTemperature

    def getNodalNormalHeatFlux(self):
        return self.nodalNormalHeatFlux

    def getNodalHeatFluxes(self):
        return (self.nodalHeatFlux_X, self.nodalHeatFlux_Y, self.nodalHeatFlux_Z)

    def applyNodalDisplacements(self, dx, dy, dz, dx_nM1, dy_nM1, dz_nM1, haloNodesDisplacements, dt):
        return

    def applyNodalHeatFluxes(self, HF_X, HF_Y, HF_Z, dt):
        return

    def applyNodalTemperatures(self, Temperature, dt):
        return

    def update(self, dt):
        return

    def steadyUpdate(self):
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

    def boundaryConditionsUpdate(self):
        return

    def exit(self):
        return

# ----------------------------------------------------------------------
#  Generic solid adjoint solver class
# ----------------------------------------------------------------------

class SolidAdjointSolver(SolidSolver):
    def __init__(self):

        SolidSolver.__init__(self)

        self.nodalAdjDisp_X = np.zeros((self.nPhysicalNodes))
        self.nodalAdjDisp_Y = np.zeros((self.nPhysicalNodes))
        self.nodalAdjDisp_Z = np.zeros((self.nPhysicalNodes))

        self.nodalAdjLoad_X = np.zeros((self.nPhysicalNodes))
        self.nodalAdjLoad_Y = np.zeros((self.nPhysicalNodes))
        self.nodalAdjLoad_Z = np.zeros((self.nPhysicalNodes))

    def applyNodalAdjointDisplacement(self, disp_adj_X, disp_adj_Y, disp_adj_Z, haloNodesDisplacements, dt):
        return
    
    def getNodalAdjointLoads(self):
        return (self.nodalAdjLoad_X, self.nodalAdjLoad_Y, self.nodalAdjDisp_Z)

# ----------------------------------------------------------------------
#  Generic fluid adjoint solver class
# ----------------------------------------------------------------------

class FluidAdjointSolver(FluidSolver):
    def __init__(self):
        
        FluidSolver.__init__(self)
        self.nodalAdjDisp_X = np.zeros((self.nPhysicalNodes))
        self.nodalAdjDisp_Y = np.zeros((self.nPhysicalNodes))
        self.nodalAdjDisp_Z = np.zeros((self.nPhysicalNodes))

        self.nodalAdjLoad_X = np.zeros((self.nPhysicalNodes))
        self.nodalAdjLoad_Y = np.zeros((self.nPhysicalNodes))
        self.nodalAdjLoad_Z = np.zeros((self.nPhysicalNodes))

    def applyNodalAdjointLoads(self, load_adj_X, load_adj_Y, load_adj_Z, haloNodesLoads, dt):
        return
    
    def getNodalAdjointDisplacement(self):
        return (self.nodalAdjDisp_X, self.nodalAdjDisp_Y, self.nodalAdjDisp_Z)