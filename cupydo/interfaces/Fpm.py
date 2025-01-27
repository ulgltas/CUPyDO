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
Fpm.py
Python interface between Fpm and CUPyDO.
Author A. Dechamps (many routines were copy-pasted from Flow.py)
"""

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import sys
import numpy as np
from cupydo.genericSolvers import FluidSolver

# ----------------------------------------------------------------------
#  Fpm Solver class
# ----------------------------------------------------------------------

class Fpm(FluidSolver):
    def __init__(self, p):
        # load the python module and initialize the solver
        module = __import__(p['cfdFile'])
        fpmP = module.getParams()
        self.__initFpm(fpmP)

        # count fsi nodes and get their positions
        self.nNodes = self.boundary.nodes.size()
        self.nHaloNode = 0
        self.nPhysicalNodes = self.nNodes - self.nHaloNode
        self.nodalInitPosX, self.nodalInitPosY, self.nodalInitPosZ = self.getNodalInitialPositions()

        # init save frequency (fsi)
        if 'saveFreq' in fpmP:
            self.saveFreq = fpmP['saveFreq']
        else:
            self.saveFreq = sys.maxsize

        # generic init
        FluidSolver.__init__(self, p)

    def __initFpm(self, p):
        """Initilize fpmw classes
        """
        import fpmw
        import tbox
        import tbox.gmsh as gmsh

        # basic checks
        if 'AoS' in p and p['AoS'] != 0 and 'Symmetry' in p:
            raise RuntimeError('Symmetry boundary condition cannot be used with nonzero angle of sideslip (AoS)!\n')
        # basic config
        if p['Format'] == 'vtk':
            try:   
               import tboxVtk
               Writer = tboxVtk.VtkExport
               print ("Found VTK libraries! Results will be saved in VTK format.\n")
            except:
               Writer = tbox.GmshExport
               print ("VTK libraries not found! Results will be saved in gmsh format.\n")
        else:
            Writer = tbox.GmshExport
            print ("Results will be saved in gmsh format.\n")

        # mesh the geometry
        self.msh = gmsh.MeshLoader(p['File'],__file__).execute(**p['Pars'])
        self.mshWriter = Writer(self.msh)

        # initialize the problem
        p['AoA'] = p['AoA']*np.pi/180 # convert to radians
        if 'AoS' in p:
            p['AoS'] = p['AoS']*np.pi/180
        else:
            p['AoS'] = 0.
        self.dynP = p['P_dyn']
        self.pbl = fpmw.Problem(self.msh, p['AoA'], p['AoS'], p['M_inf'], p['S_ref'], p['c_ref'], p['x_ref'], p['y_ref'], p['z_ref'], bool(p['sym']))    
        
        # define bodies and identify f/s boundary
        for i in range(len(p['Wings'])):
            bnd = fpmw.Body(self.msh, p['Wings'][i], [p['Tes'][i]], p['LengthWake'][i])
            self.pbl.add(bnd)
            if p['Wings'][i] == p['Fsi']:
                self.boundary = bnd
        
        # initialize the AIC builder
        self.aic = fpmw.Builder(self.pbl)

        # initialize the solver
        self.solver = fpmw.Solver(self.aic)    
        
    def run(self, t1, t2):
        """Run the solver for one steady (time) iteration
        """
        self.solver.run()
        self.__setCurrentState()
        return True
    
    def __setCurrentState(self):
        """Compute nodal forces from nodal normalized forces
        """
        i = 0
        for n in self.boundary.nodes:
            self.nodalLoad_X[i] = -self.dynP * self.boundary.cLoadX[i]
            self.nodalLoad_Y[i] = -self.dynP * self.boundary.cLoadY[i]
            self.nodalLoad_Z[i] = -self.dynP * self.boundary.cLoadZ[i]
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

    def applyNodalDisplacements(self, dx, dy, dz, dt, haloNodesDisplacements):
        """Apply displacements coming from solid solver to f/s interface after saving
        """
        for i in range(self.boundary.nodes.size()):
            self.boundary.nodes[i].pos[0] = self.nodalInitPosX[i] + dx[i]
            self.boundary.nodes[i].pos[1] = self.nodalInitPosY[i] + dy[i]
            self.boundary.nodes[i].pos[2] = self.nodalInitPosZ[i] + dz[i]

    def meshUpdate(self, nt):
        """Deform the mesh using linear elasticity equations
        """
        self.aic.run()
        
    def save(self, nt):
        """Save data on disk at each converged timestep
        """
        self.solver.save(self.mshWriter, nt)
        self.mshWriter.save(self.msh.name + "_" + str(nt))

    def initRealTimeData(self):
        """Initialize history file
        """
        histFile = open('FlowHistory.dat', 'w')
        histFile.write('{0:>12s}   {1:>12s}   {2:>12s}   {3:>12s}   {4:>12s}\n'.format('Time', 'FSI_Iter', 'C_Lift', 'C_Drag', 'C_Moment'))
        histFile.close()

    def saveRealTimeData(self, time, nFSIIter):
        """Save data at each fsi iteration
        """
        # history at each iteration
        histFile = open('FlowHistory.dat', 'a')
        histFile.write('{0:12.6f}   {1:12d}   {2:12.6f}   {3:12.6f}   {4:12.6f}\n'.format(time, nFSIIter, self.solver.Cl, self.solver.Cd, self.solver.Cm))
        histFile.close()
        # full solution at user-defined frequency
        if np.mod(nFSIIter+1, self.saveFreq) == 0:
            self.solver.save(self.mshWriter, 1000000+int(nFSIIter+1)//int(self.saveFreq))

    def printRealTimeData(self, time, nFSIIter):
        """Print data on screen at the end of fsi simulation
        """
        print ('[Flow lift, drag, moment]: {0:6.3f}, {1:6.4f}, {2:6.3f}'.format(self.solver.Cl, self.solver.Cd, self.solver.Cm))
        print ('')
    
    def exit(self):
        """
        Exit the Flow solver
        """
        del self.solver
        del self.aic
        del self.pbl
        del self.mshWriter
        del self.msh
