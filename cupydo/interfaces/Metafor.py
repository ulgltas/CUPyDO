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

# %% Nodal Load class

class NLoad(object):
    """
    Class representing the nodal loads/temperatures
    """
    def __init__(self,val1,t1,val2,t2):

        self.t1 = t1
        self.t2 = t2
        self.val1 = val1
        self.val2 = val2

    def __call__(self,time):

        return self.val1+(time-self.t1)/(self.t2-self.t1)*(self.val2-self.val1)

    def nextstep(self):

        self.t1 = self.t2
        self.val1 = self.val2

# %% Interface Between Metafor and CUPyDO
               
class Metafor(SolidSolver):

    def __init__(self,param):

        titlePrint('Initializing Metafor')

        input = dict()
        module = importlib.import_module(param['csdFile'])
        self.metafor = module.getMetafor(input)
        domain = self.metafor.getDomain()
        input = module.params(input)

        # Defines some internal variables

        self.Fnods = dict()
        self.Tnods = dict()
        self.neverRun = True
        self.reload = True

        # Defines some internal variables

        self.exporter = input['exporter']
        loadingset = domain.getLoadingSet()
        self.computationType = param['compType']
        self.tsm = self.metafor.getTimeStepManager()
        self.FSI = domain.getGeometry().getGroupSet()(input['bndno'])
        self.nNodes = self.FSI.getNumberOfMeshPoints()
        self.nPhysicalNodes = self.nNodes

        # Creates the nodal load container

        for i in range(self.nNodes):
            
            load = list()
            temp = NLoad(0,0,0,0)
            node = self.FSI.getMeshPoint(i)
            self.Fnods[node.getNo()] = load
            self.Tnods[node.getNo()] = temp

            fct = w.PythonOneParameterFunction(temp)
            # loadingset.define(node,w.Field1D(w.TO,w.AB),1,fct)
            # loadingset.define(node,w.Field1D(w.TO,w.RE),1,fct,0)

            for F in [w.TX,w.TY,w.TZ]:

                load.append(NLoad(0,0,0,0))
                fct = w.PythonOneParameterFunction(load[-1])
                loadingset.define(node,w.Field1D(F,w.GF1),1,fct)

        # Creates the array for external communication

        self.vel = np.zeros((self.nNodes,3))
        self.acc = np.zeros((self.nNodes,3))
        self.dis = np.zeros((self.nNodes,3))

        # Initialization of domain and output

        SolidSolver.__init__(self)
        self.metafor.getDomain().build()
        self.__setCurrentState(True)
        self.save()

        # Manages time step restart functions

        self.mfac = w.MemoryFac()
        self.metaFac = w.MetaFac(self.metafor)
        self.metaFac.mode(False,False,True)
        self.metaFac.save(self.mfac)
        self.tsm.setVerbose(False)

# %% Calculates One Time Step

    def run(self,t1,t2):
        """
        Computes a time increment and/or load previous state
        """

        if(self.neverRun):

            self.tsm.setInitialTime(t1,t2-t1)
            self.tsm.setNextTime(t2,0,t2-t1)
            ok = self.metafor.getTimeIntegration().integration()
            self.neverRun = False

        else:

            if self.reload: self.tsm.removeLastStage()
            self.tsm.setNextTime(t2,0,t2-t1)
            ok = self.metafor.getTimeIntegration().restart(self.mfac)

        self.__setCurrentState(False)
        self.reload = True
        return ok

# %% Gets Nodal Values

    def __setCurrentState(self,initialize):
        """
        Save the current nodal states in vectors
        """

        for i in range(self.nNodes):
            node = self.FSI.getMeshPoint(i)

            self.nodalVel_X[i] = node.getValue(w.Field1D(w.TX,w.GV))
            self.nodalVel_Y[i] = node.getValue(w.Field1D(w.TY,w.GV))
            self.nodalVel_Z[i] = node.getValue(w.Field1D(w.TZ,w.GV))

            self.nodalDisp_X[i] = node.getValue(w.Field1D(w.TX,w.RE))
            self.nodalDisp_Y[i] = node.getValue(w.Field1D(w.TY,w.RE))
            self.nodalDisp_Z[i] = node.getValue(w.Field1D(w.TZ,w.RE))

            if not initialize:

                HF_X_extractor = w.IFNodalValueExtractor(node,w.IF_FLUX_X)
                HF_Y_extractor = w.IFNodalValueExtractor(node,w.IF_FLUX_Y)
                HF_Z_extractor = w.IFNodalValueExtractor(node,w.IF_FLUX_Z)

                self.nodalHeatFlux_X[i] = HF_X_extractor.extract()[0]
                self.nodalHeatFlux_Y[i] = HF_Y_extractor.extract()[0]
                self.nodalHeatFlux_Z[i] = HF_Z_extractor.extract()[0]

    def getNodalInitialPositions(self):
        """
        Return the position vector of the solver interface
        """

        pos = np.zeros((3,self.nNodes))
        for i in range(self.nNodes):

            node = self.FSI.getMeshPoint(i)
            pos[0,i] = node.getValue(w.Field1D(w.TX,w.AB))
            pos[1,i] = node.getValue(w.Field1D(w.TY,w.AB))
            pos[2,i] = node.getValue(w.Field1D(w.TZ,w.AB))

        return pos

    def getNodalIndex(self,index):
        """
        Return the index of the index-th interface node in the load vector
        """

        node = self.FSI.getMeshPoint(index)
        return node.getNo()

# %% Set Nodal Loads

    def applyNodalLoads(self,loadX,loadY,loadZ,time,*_):
        """
        Apply the load boundary conditions on the mesh
        """

        for i in range(self.nNodes):

            node = self.FSI.getMeshPoint(i)
            fx,fy,fz = self.Fnods[node.getNo()]
            fx.val2 = loadX[i]
            fy.val2 = loadY[i]
            fz.val2 = loadZ[i]
            fx.t2 = time
            fy.t2 = time
            fz.t2 = time

    def applyNodalTemperatures(self,temp,time):
        """
        Apply the temperature boundary conditions on the mesh
        """

        for i in range(self.nNodes):

            node = self.FSI.getMeshPoint(i)
            temp = self.Tnods[node.getNo()]
            temp.val2 = temp[i]
            temp.t2 = time

# %% Other Functions

    def update(self):
        """
        Save the current state in the RAM and update the load
        """

        for F in self.Fnods.values(): [F[i].nextstep() for i in range(3)]
        for T in self.Tnods.values(): T.nextstep()
        self.metaFac.save(self.mfac)
        SolidSolver.update(self)
        self.reload = False

    def steadyUpdate(self):
        """
        Save the current state in the RAM and update the load
        """

        self.update()

    def save(self):
        """
        Save the curent state to the disk
        """

        if self.exporter is not None: self.exporter.execute()
    
    def exit(self):
        """
        Exit the solid solver
        """

        titlePrint('Exit Metafor')