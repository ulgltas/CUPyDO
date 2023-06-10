#!/usr/bin/env python3
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

Modal.py
Python interface between a modal solver and CUPyDO.
Adrien Crovato, Mariano Sanchez Martinez

'''

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import numpy as np
from ..utilities import titlePrint
from ..genericSolvers import SolidSolver

# ----------------------------------------------------------------------
#  Modal solver interface class
# ----------------------------------------------------------------------

class Modal(SolidSolver):
    """
    Modal interface for CUPyDO
    """

    def __init__(self, _module, _computationType):

        titlePrint('Initialize Modal Interface')
        # load the python module and initialize modal solver
        module = __import__(_module)
        self.initModal(module.getParams())

        # get number of nodes
        self.nNodes = self.solver.nNodes
        self.nHaloNodes = 0
        self.nPhysicalNodes = self.nNodes - self.nHaloNodes

        # initialize
        SolidSolver.__init__(self)
        self.computationType = _computationType
        self.setInitialDisplacements()
        self.initRealTimeData()

    def initModal(self, p):
        """Initialize ModalSolver classes
        Adrien Crovato
        """
        import modali
        # initialize solver
        self.solver = modali.modali(p['nm'])
        self.solver.setMatrices(p['M_q'], p['C_q'], p['K_q']) # modal matrices
        self.solver.readModes(p['File']) # config file
        self.solver.setInitial(p['x_i'], p['v_i'], p['f_i']) # initial conditions
        self.solver.setExtractor(p['Extractors']) # extractor list

    def setInitialDisplacements(self):
        """Set initial displacements
        Adrien Crovato
        """
        self.__setCurrentState()

    def run(self, t1, t2):
        """Run the solver between t1 and t2
        Adrien Crovato
        """
        if self.computationType == 'steady':
            self.solver.runStatic()
        else:
            self.solver.runDynamic(t1, t2)

        self.__setCurrentState()

    def __setCurrentState(self):
        """Update displacements
        Adrien Crovato
        """
        self.nodalDisp_X = self.solver.dispX
        self.nodalDisp_Y = self.solver.dispY
        self.nodalDisp_Z = self.solver.dispZ
           
    def applyNodalLoads(self, load_X, load_Y, load_Z, time, haloNodesLoads = {}):
        """Update the loads
        Adrien Crovato
        """
        self.solver.updateLoads(load_X[0:self.nPhysicalNodes], load_Y[0:self.nPhysicalNodes], load_Z[0:self.nPhysicalNodes])            
            
    def getNodalInitialPositions(self):
        """Return initial nodal positions
        Adrien Crovato, Mariano Sanchez Martinez
        """
        # Numpy arrays have to be C-contiguous in order to be passed through MPI
        return (np.ascontiguousarray(self.solver.nodalCoord_X),
            np.ascontiguousarray(self.solver.nodalCoord_Y),
            np.ascontiguousarray(self.solver.nodalCoord_Z))

    def getNodalIndex(self, iVertex):
        """Get index of nodes
        Adrien Crovato
        """
        return self.solver.nodalGlobalIndex[iVertex]

    def initRealTimeData(self):
        """Initialize results files (history and extractors)
        Adrien Crovato
        """
        # History
        histFile = open('ModalHistory.dat', "w")
        histFile.write('{0:>12s}   {1:>12s}'.format("Time", "FSI_Iter"))
        for i in range(0, self.solver.nModes):
            histFile.write('   {0:>12s}'.format('y_'+str(i)))
        for i in range(0, self.solver.nModes):
            histFile.write('   {0:>12s}'.format('fq_'+str(i)))
        histFile.write('\n')
        histFile.close()
        # Nodal displacements
        solFile = open('NodalDisplacement.dat', "w")
        solFile.write('{0:>12s}   {1:>12s}'.format("Time", "FSI_Iter"))
        for gidx in self.solver.extractor:
            solFile.write('   {0:>12s}   {1:>12s}   {2:>12s}'.format('x_'+str(gidx), 'y_'+str(gidx), 'z_'+str(gidx)))
        solFile.write('\n')
        solFile.close()


    def saveRealTimeData(self, time, nFSIIter):
        """Write results to file (history and extractor)
        Adrien Crovato
        """
        # History
        histFile = open('ModalHistory.dat', "a")
        histFile.write("{0:12.6f}   {1:12d}".format(time, nFSIIter))
        for i in range(0, self.solver.nModes):
            histFile.write('   {0:12.6f}'.format(self.solver.y0[i]))
        for i in range(0, self.solver.nModes):
            histFile.write('   {0:12.6f}'.format(self.solver.fq[i]))
        histFile.write('\n')
        histFile.close()
        # Nodal displacements
        solFile = open('NodalDisplacement.dat', "a")
        solFile.write("{0:12.6f}   {1:12d}".format(time, nFSIIter))
        for lidx in list(self.solver.extractor.values()):
            solFile.write('   {0:12.6f}   {1:12.6f}   {2:12.6f}'.format(self.nodalDisp_X[lidx], self.nodalDisp_Y[lidx], self.nodalDisp_Z[lidx]))
        solFile.write('\n')
        solFile.close()

    def save(self):
        """Save data on disk at each converged timestep
        Adrien Crovato
        """
        self.solver.write('deformed_modes')

    def exit(self):
        """Say bye
        Adrien Crovato.
        """
        del self.solver
        
        titlePrint("Exit Modal Interface")
