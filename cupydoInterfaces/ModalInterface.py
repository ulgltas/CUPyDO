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
        self.nNodes = self.modal.solver.nNodes
        self.nHaloNodes = 0
        self.nPhysicalNodes = self.nNodes - self.nHaloNodes

        # initialize
        SolidSolver.__init__(self)
        self.computationType = _computationType
        self.__setCurrentState()
        self.initRealTimeData()

    def setInitialDisplacements(self):
        """
        Des.
        """
        self.__setCurrentState()

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
        self.modal.solver.updateLoads(load_X[0:self.nPhysicalNodes], load_Y[0:self.nPhysicalNodes], load_Z[0:self.nPhysicalNodes])            
            
    def getNodalInitialPositions(self):
        """
        Des.
        """
        # Numpy arrays have to be C-contiguous in order to be passed through MPI
        return (np.ascontiguousarray(self.modal.solver.nodalCoord_X),
            np.ascontiguousarray(self.modal.solver.nodalCoord_Y),
            np.ascontiguousarray(self.modal.solver.nodalCoord_Z))

    def getNodalIndex(self, iVertex):
        """
        Des.
        """
        return self.modal.solver.nodalGlobalIndex[iVertex]

    def initRealTimeData(self):
        """Initialize results files
        """
        # History
        histFile = open('ModalHistory.dat', "w")
        histFile.write('{0:>12s}   {1:>12s}'.format("Time", "FSI_Iter"))
        for i in range(0, self.modal.solver.nModes):
            histFile.write('   {0:>12s}'.format('y_'+str(i)))
        for i in range(0, self.modal.solver.nModes):
            histFile.write('   {0:>12s}'.format('fq_'+str(i)))
        histFile.write('\n')
        histFile.close()
        # Nodal displacements
        solFile = open('NodalDisplacement.dat', "w")
        solFile.write('{0:>12s}   {1:>12s}'.format("Time", "FSI_Iter"))
        for gidx in self.modal.solver.extractor:
            solFile.write('   {0:>12s}   {1:>12s}   {2:>12s}'.format('x_'+str(gidx), 'y_'+str(gidx), 'z_'+str(gidx)))
        solFile.write('\n')
        solFile.close()


    def saveRealTimeData(self, time, nFSIIter):
        """Write results to file
        """
        # History
        histFile = open('ModalHistory.dat', "a")
        histFile.write("{0:12.6f}   {1:12d}".format(time, nFSIIter))
        for i in range(0, self.modal.solver.nModes):
            histFile.write('   {0:12.6f}'.format(self.modal.solver.y0[i]))
        for i in range(0, self.modal.solver.nModes):
            histFile.write('   {0:12.6f}'.format(self.modal.solver.fq[i]))
        histFile.write('\n')
        histFile.close()
        # Nodal displacements
        solFile = open('NodalDisplacement.dat', "a")
        solFile.write("{0:12.6f}   {1:12d}".format(time, nFSIIter))
        for lidx in self.modal.solver.extractor.values():
            solFile.write('   {0:12.6f}   {1:12.6f}   {2:12.6f}'.format(self.nodalDisp_X[lidx], self.nodalDisp_Y[lidx], self.nodalDisp_Z[lidx]))
        solFile.write('\n')
        solFile.close()

    def save(self):
        """
        Save data on disk at each converged timestep
        """
        self.modal.solver.write('deformed_modes')

    def exit(self):
        """
        Des.
        """
        print("***************************** Exit modal interface *****************************")
