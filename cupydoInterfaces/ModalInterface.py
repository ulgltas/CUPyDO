#!/usr/bin/env python
# -*- coding: latin-1; -*-

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

ModalInterface.py
Python interface between a modal solver and CUPyDO.
Huseyin Guner, Adrien Crovato 

'''

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import numpy as np
from cupydo.genericSolvers import SolidSolver

# ----------------------------------------------------------------------
#  Modal solver interface class
# ----------------------------------------------------------------------

class ModalInterface(SolidSolver):
    """
    Modal interface for CUPyDO
    """

    def __init__(self, _module, _computationType):
        print("\n***************************** Initialize modal interface *****************************\n")
        # load the python module
        module = __import__(_module)
        self.modal = module.getModal()

        # get number of nodes
        self.nNodes = self.modal.nNodes
        self.nHaloNodes = 0
        self.nPhysicalNodes = self.nNodes - self.nHaloNodes

        # initialize
        SolidSolver.__init__(self)
        self.computationType = _computationType
        self.__setCurrentState()
        self.initRealTimeData()

    def run(self, t1, t2):
        """
        Des.
        """
        if self.computationType == 'steady':
            self.modal.solver.runStatic()
        else:
            self.modal.solver.runDynamic(t1, t2)

        self.__setCurrentState()

    def __setCurrentState(self):
        """
        Des.
        """
        self.nodalDisp_X = self.modal.solver.dispX
        self.nodalDisp_Y = self.modal.solver.dispY
        self.nodalDisp_Z = self.modal.solver.dispZ          
           
    def applyNodalLoads(self, load_X, load_Y, load_Z, time):
        """
        """
        self.modal.updateLoads(load_X, load_Y, load_Z)            
            
    def getNodalInitialPositions(self):
        """
        Des.
        """
        return (self.modal.nodalCoord_X, self.modal.nodalCoord_Y, self.modal.nodalCoord_Z)

    def getNodalIndex(self, iVertex):
        """
        Des.
        """
        return self.modal.nodalGlobalIndex[iVertex]

    #def update(self):
    #    """
    #    Des.
    #    """
    #    SolidSolver.update(self)

    def initRealTimeData(self):
        """
        """
        histFile = open('ModalHistory.dat', "w")
        histFile.write('{0:>12s}   {1:>12s}'.format("Time", "FSI_Iter"))
        for i in range(0, self.modal.solver.nModes)
            histFile.write('   {0:>12s}'.format('y_'+str(i)))
        for i in range(0, self.modal.solver.nModes)
            histFile.write('   {0:>12s}'.format('fq_'+str(i)))
        histFile.write('\n')
        histFile.close()

    def saveRealTimeData():
        """
        """
        histFile = open('ModalHistory.dat', "a")
        histFile.write('{0:12.6f}   {1:12d}'.format("Time", "FSI_Iter"))
        for i in range(0, self.modal.solver.nModes)
            histFile.write('   {0:12.6f}'.format(self.modal.solver.y0[i]))
        for i in range(0, self.modal.solver.nModes)
            histFile.write('   {0:12.6f}'.format(self.modal.solver.fq[i]))
        histFile.write('\n')
        histFile.close()

    def exit(self):
        """
        Des.
        """
        print("***************************** Exit modal interface *****************************")
