#!/usr/bin/env python
# -*- coding: utf-8; -*-

"""
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

Dart.py
Python interface between DART and CUPyDO.
Authors A. Crovato
"""

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import sys
import numpy as np
from ..genericSolvers import FluidSolver

# ----------------------------------------------------------------------
#  DartSolver class
# ----------------------------------------------------------------------

class Dart(FluidSolver):
    def __init__(self, _module, _nthreads):
        # load the python module and initialize the solver
        module = __import__(_module)
        params = module.getParams()
        params['Threads'] = _nthreads
        from dart.api.core import initDart
        _, self.qinf, self.msh, self.writer, self.morpher, _, self.boundary, self.solver, _ = initDart(params, scenario='aerostructural', task='analysis')

        # count fsi nodes and get their positions
        self.nNodes = self.boundary.nodes.size()
        self.nHaloNode = 0
        self.nPhysicalNodes = self.nNodes - self.nHaloNode
        self.nodalInitPosX, self.nodalInitPosY, self.nodalInitPosZ = self.getNodalInitialPositions()

        # init save frequency (fsi)
        if 'saveFreq' in params:
            self.saveFreq = params['saveFreq']
        else:
            self.saveFreq = sys.maxsize

        # generic init
        FluidSolver.__init__(self)
        
    def run(self, t1, t2):
        """Run the solver for one steady (time) iteration.
        """
        status = self.solver.run()
        if status > 1:
            raise RuntimeError('DART solver diverged!\n')
        self.__setCurrentState()
    
    def __setCurrentState(self):
        """Compute nodal forces from nodal normalized forces
        """
        i = 0
        for n in self.boundary.nodes:
            self.nodalLoad_X[i] = self.qinf * self.boundary.nLoads[i][0]
            self.nodalLoad_Y[i] = self.qinf * self.boundary.nLoads[i][1]
            self.nodalLoad_Z[i] = self.qinf * self.boundary.nLoads[i][2]
            i += 1

    def getNodalInitialPositions(self):
        """Get the initial position of each node
        """
        x0 = np.zeros(self.nPhysicalNodes)
        y0 = np.zeros(self.nPhysicalNodes)
        z0 = np.zeros(self.nPhysicalNodes)
        for i in range(self.boundary.nodes.size()):
            n = self.boundary.nodes[i]               
            x0[i] = n.pos[0]
            y0[i] = n.pos[1]
            z0[i] = n.pos[2]

        return (x0, y0, z0)

    def getNodalIndex(self, iVertex):
        """Get index of each node
        """
        no = self.boundary.nodes[iVertex].no
        return no

    def applyNodalDisplacements(self, dx, dy, dz, dx_nM1, dy_nM1, dz_nM1, haloNodesDisplacements, time):
        """Apply displacements coming from solid solver to f/s interface after saving
        """
        self.morpher.savePos()
        for i in range(self.boundary.nodes.size()):
            self.boundary.nodes[i].pos[0] = self.nodalInitPosX[i] + dx[i]
            self.boundary.nodes[i].pos[1] = self.nodalInitPosY[i] + dy[i]
            self.boundary.nodes[i].pos[2] = self.nodalInitPosZ[i] + dz[i]

    def meshUpdate(self, nt):
        """Deform the mesh using linear elasticity equations
        """
        self.morpher.deform()
        
    def save(self, nt):
        """Save data on disk at each converged timestep
        """
        self.solver.save(self.writer, nt)
        self.writer.save(self.msh.name + "_" + str(nt))

    def initRealTimeData(self):
        """Initialize history file
        """
        histFile = open('DartHistory.dat', 'w')
        histFile.write('{0:>12s}   {1:>12s}   {2:>12s}   {3:>12s}   {4:>12s}\n'.format('Time', 'FSI_Iter', 'C_Lift', 'C_Drag', 'C_Moment'))
        histFile.close()

    def saveRealTimeData(self, time, nFSIIter):
        """Save data at each fsi iteration
        """
        # history at each iteration
        histFile = open('DartHistory.dat', 'a')
        histFile.write('{0:12.6f}   {1:12d}   {2:12.6f}   {3:12.6f}   {4:12.6f}\n'.format(time, nFSIIter, self.solver.Cl, self.solver.Cd, self.solver.Cm))
        histFile.close()
        # full solution at user-defined frequency
        if np.mod(nFSIIter+1, self.saveFreq) == 0:
            self.solver.save(self.writer, 1000000+int(nFSIIter+1)//int(self.saveFreq))

    def printRealTimeData(self, time, nFSIIter):
        """Print data on screen at the end of fsi simulation
        """
        print('[DART lift, drag, moment]: {0:6.3f}, {1:6.4f}, {2:6.3f}'.format(self.solver.Cl, self.solver.Cd, self.solver.Cm))
        print('')
    
    def exit(self):
        """Clear memory
        """
        del self.boundary
        del self.solver
        del self.morpher
        del self.writer
        del self.msh
