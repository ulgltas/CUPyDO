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

Metafor.py
Python interface between the wrapper of Metafor and CUPyDO.
Authors R. BOMAN, M.L. CERQUAGLIA, D. THOMAS

'''

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

from ..utilities import titlePrint
from ..genericSolvers import SolidSolver
import numpy as np
import wrap as w
import importlib

# ----------------------------------------------------------------------
#  Nodal Load class
# ----------------------------------------------------------------------

class NLoad(object):
    """
    Class representing the nodal forces
    """
    def __init__(self, val):
        self.val = float(val)

    def __call__(self, time):
        return float(self.val)

# ----------------------------------------------------------------------
#  Metafor solver interface class
# ----------------------------------------------------------------------
               
class Metafor(SolidSolver):

    def __init__(self, p):

        titlePrint('Initializing Metafor')

        parm = dict()
        module = importlib.import_module(p['csdFile'])
        self.metafor = module.getMetafor(parm)
        domain = self.metafor.getDomain()

        # Defines some internal variables

        self.reload = True
        self.__dict__.update(parm)

        self.regime = p['regime']
        self.interpType = p['interpType']

        self.thermal = p['thermal']
        self.mechanical = p['mechanical']

        # Get Metafor Python objects

        geometry = domain.getGeometry()
        loadingset = domain.getLoadingSet()
        self.tsm = self.metafor.getTimeStepManager()

        self.dim = geometry.getDimension().getNdim()
        self.nNodes = self.FSI.getNumberOfMeshPoints()
        self.nPhysicalNodes = self.nNodes

        # Creates the conservative nodal load container

        if self.interpType == 'conservative':

            self.Fnods = dict()
            for i in range(self.nNodes):
                
                load = list()
                node = self.FSI.getMeshPoint(i)
                self.Fnods[node.getNo()] = load

                for F in w.TX, w.TY, w.TZ:

                    load.append(NLoad(0))
                    fct = w.PythonOneParameterFunction(load[-1])
                    loadingset.define(node,w.Field1D(F, w.GF1), 1, fct)

        # Initialization of domain and output

        SolidSolver.__init__(self)

        self.metafor.getTimeStepManager().setInitialTime(0, 1)
        self.metafor.getTimeIntegration().initialise()
        self.__setCurrentState()
        self.save()

        # Manages time step restart functions

        self.mfac = w.MemoryFac()
        self.metaFac = w.MetaFac(self.metafor)
        self.metaFac.mode(False,False, True)
        self.metaFac.save(self.mfac)
        self.tsm.setVerbose(False)

    def run(self, t1, t2):
        """
        Computes a time increment and/or load previous state
        """

        if self.reload: self._wayBack()
        self.tsm.setNextTime(t2, 0, t2-t1)

        ok = self.metafor.getTimeIntegration().restart(self.mfac)

        if ok: self.__setCurrentState()
        self.reload = True
        return ok

    def __setCurrentState(self):
        """
        Save the current nodal states in vectors
        """

        for i in range(self.nNodes):
            node = self.FSI.getMeshPoint(i)

            if self.mechanical:

                self.nodalVel_X[i, 0] = node.getValue(w.Field1D(w.TX, w.GV))
                self.nodalVel_Y[i, 0] = node.getValue(w.Field1D(w.TY, w.GV))
                self.nodalVel_Z[i, 0] = node.getValue(w.Field1D(w.TZ, w.GV))

                self.nodalDisp_X[i, 0] = node.getValue(w.Field1D(w.TX, w.RE))
                self.nodalDisp_Y[i, 0] = node.getValue(w.Field1D(w.TY, w.RE))
                self.nodalDisp_Z[i, 0] = node.getValue(w.Field1D(w.TZ, w.RE))

            if self.thermal:
                self.nodalTemperature[i] = node.getValue(w.Field1D(w.TO, w.AB)) + node.getValue(w.Field1D(w.TO, w.RE))

    def getNodalInitialPositions(self):
        """
        Return the position vector of the solver interface
        """

        pos = np.zeros((3, self.nNodes))
        for i in range(self.nNodes):

            node = self.FSI.getMeshPoint(i)
            pos[0,i] = node.getValue(w.Field1D(w.TX, w.AB))
            pos[1,i] = node.getValue(w.Field1D(w.TY, w.AB))
            pos[2,i] = node.getValue(w.Field1D(w.TZ, w.AB))

        return pos

    def getNodalIndex(self,index):
        """
        Return the index of the index-th interface node in the load vector
        """

        node = self.FSI.getMeshPoint(index)
        return node.getNo()

    def applyNodalForce(self, load_X, load_Y, load_Z, dt, haloNodesLoads):
        """
        Apply the conservative load boundary conditions
        """

        result = np.transpose([load_X, load_Y, load_Z])[:self.nNodes]
        for i in range(self.nNodes):

            node = self.FSI.getMeshPoint(i)
            fx,fy,fz = self.Fnods[node.getNo()]
            fx.val = result[i][0][0]
            fy.val = result[i][0][1]
            fz.val = result[i][0][2]


    def applyNodalStress(self, load_XX, load_YY, load_ZZ, load_XY, load_XZ, load_YZ, dt, haloNodesLoads):
        """
        Apply the consistent load boundary conditions
        """

        result = np.transpose([load_XX, load_YY, load_ZZ, load_XY, load_XZ, load_YZ])[:self.nNodes]
        for i in range(self.nNodes):

            node = self.FSI.getMeshPoint(i)
            self.interactionM.setNodTensor3D(node, *result[i])

    def applyNodalHeatFluxes(self, HF_X, HF_Y, HF_Z, dt, haloNodesHeatFlux):
        """
        Apply the consistent heat flux boundary conditions
        """

        result = np.transpose([HF_X, HF_Y, HF_Z])[:self.nNodes]
        for i in range(self.nNodes):

            node = self.FSI.getMeshPoint(i)
            self.interactionT.setNodVector(node, *result[i])

    def update(self):
        """
        Save the current state in the RAM and update the load
        """

        self.metaFac.save(self.mfac)
        SolidSolver.update(self)
        self.reload = False

    def _wayBack(self):
        '''
        Revert back the solver to its last converged FSI state
        '''

        # Restart if either the step or the time has advanced

        if (self.metafor.getCurrentStepNo() > self.mfac.getStepNo()) or \
           (self.metafor.getLastTime() >= self.metafor.getCurrentTime()):

            self.metafor.load(self.mfac)

        # Remove the last stage if more than one stage has been stored

        if not (self.metafor.getStageManager().getCurNumStage() < 0) \
           and (self.metafor.getStageManager().getNumbOfStage() > 1):

            self.metafor.getTimeStepManager().removeLastStage()

    def steadyUpdate(self):
        """
        Save the current state in the RAM and update the load
        """

        self.metaFac.save(self.mfac)
        self.reload = True

    def save(self):
        """
        Save the curent state to the disk
        """

        if self.exporter is not None:
            self.exporter.write()
    
    def exit(self):
        """
        Exit the solid solver
        """

        titlePrint('Exit Metafor')
