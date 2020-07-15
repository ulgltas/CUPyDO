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

Flow.py
Python interface between Flow and CUPyDO.
Authors A. Crovato
"""
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
from builtins import str
from builtins import range
from past.utils import old_div
import sys
import numpy as np
from ..genericSolvers import FluidSolver

# ----------------------------------------------------------------------
#  FlowSolver class
# ----------------------------------------------------------------------

class Flow(FluidSolver):
    def __init__(self, _module, _nthreads):
        # load the python module and initialize the solver
        module = __import__(_module)
        floP = module.getParams()
        self.__initFlow(floP, _nthreads)

        # count fsi nodes and get their positions
        self.nNodes = self.boundary.nodes.size()
        self.nHaloNode = 0
        self.nPhysicalNodes = self.nNodes - self.nHaloNode
        self.nodalInitPosX, self.nodalInitPosY, self.nodalInitPosZ = self.getNodalInitialPositions()

        # init save frequency (fsi)
        if 'saveFreq' in floP:
            self.saveFreq = floP['saveFreq']
        else:
            self.saveFreq = sys.maxsize

        # generic init
        FluidSolver.__init__(self)

    def __initFlow(self, p, _nthreads):
        """Initilize flow classes
        Adrien Crovato
        """
        from . import flow
        import tbox
        import tbox.gmsh as gmsh
        from tbox.solvers import LinearSolver

        # basic checks
        if p['Dim'] != 2 and p['Dim'] != 3:
            raise RuntimeError('Problem dimension should be 2 or 3, but ' + p['Dim'] + ' was given!\n')
        if p['Dim'] == 2:
            if 'AoS' in p and p['AoS'] != 0:
                raise RuntimeError('Angle of sideslip (AoS) should be zero for 2D problems!\n')
            if 'Symmetry' in p:
                raise RuntimeError('Symmetry boundary condition cannot be used for 2D problems!\n')
        else:
            if 'AoS' in p and p['AoS'] != 0 and 'Symmetry' in p:
                raise RuntimeError('Symmetry boundary condition cannot be used with nonzero angle of sideslip (AoS)!\n')
        # basic config
        if p['Format'] == 'vtk':
            try:   
               import tboxVtk
               Writer = tboxVtk.VtkExport
               print("Found VTK libraries! Results will be saved in VTK format.\n")
            except:
               Writer = tbox.GmshExport
               print("VTK libraries not found! Results will be saved in gmsh format.\n")
        else:
            Writer = tbox.GmshExport
            print("Results will be saved in gmsh format.\n")

        # mesh the geometry
        self.msh = gmsh.MeshLoader(p['File'],__file__).execute(**p['Pars'])
        if p['Dim'] == 2:
            mshCrck = tbox.MshCrack(self.msh, p['Dim'])
            mshCrck.setCrack(p['Wake'])
            mshCrck.addBoundaries([p['Fluid'], p['Farfield'][-1], p['Wing']])
            mshCrck.run()
        else:
            for i in range(0, len(p['Wakes'])):
                mshCrck = tbox.MshCrack(self.msh, p['Dim'])
                mshCrck.setCrack(p['Wakes'][i])
                mshCrck.addBoundaries([p['Fluid'], p['Farfield'][-1], p['Wings'][i]])
                if 'Fuselage' in p:
                    mshCrck.addBoundaries([p['Fuselage']])
                if 'Symmetry' in p:
                    mshCrck.addBoundaries([p['Symmetry']])
                mshCrck.setExcluded(p['WakeTips'][i])
                mshCrck.run()
        tbox.GmshExport(self.msh).save(self.msh.name)
        del mshCrck
        self.mshWriter = Writer(self.msh)

        # initialize mesh deformation handler
        self.mshDef = tbox.MshDeform(self.msh, p['Dim'])
        self.mshDef.nthreads = _nthreads
        self.mshDef.setField(p['Fluid'])
        self.mshDef.addFixed(p['Farfield'])
        if p['Dim'] == 2:
            self.mshDef.addMoving([p['Wing']])
            self.mshDef.addInternal([p['Wake'], p['Wake']+'_'])
        else:
            if 'Fuselage' in p:
                self.mshDef.addFixed(p['Fuselage'])
            self.mshDef.setSymmetry(p['Symmetry'], 1)
            for i in range(0, len(p['Wings'])):
                if p['Wings'][i] == p['Fsi']:
                    self.mshDef.addMoving([p['Wings'][i]])
                else:
                    self.mshDef.addFixed([p['Wings'][i]])
                self.mshDef.addInternal([p['Wakes'][i], p['Wakes'][i]+'_'])
        self.mshDef.initialize()

        # initialize the problem
        p['AoA'] = p['AoA']*np.pi/180 # convert to radians
        if 'AoS' in p:
            p['AoS'] = p['AoS']*np.pi/180
        else:
            p['AoS'] = 0.
        phiInfFun = flow.F0PsPhiInf(p['Dim'], p['AoA'], p['AoS'])
        velInfFun = flow.F1ElVi(p['Dim'], p['AoA'], p['AoS'])
        self.dynP = p['P_dyn']    
        pbl = flow.Problem(self.msh, p['Dim'], p['AoA'], p['AoS'], p['M_inf'], p['S_ref'], p['c_ref'], p['x_ref'], p['y_ref'], p['z_ref'])

        # add medium
        if p['M_inf'] == 0:
            pbl.set(flow.Medium(self.msh, p['Fluid'], flow.F0ElRhoL(), flow.F0ElMachL(), flow.F0ElCpL(), phiInfFun))
        else:
            pbl.set(flow.Medium(self.msh, p['Fluid'], flow.F0ElRho(p['M_inf'], p['M_crit']), flow.F0ElMach(p['M_inf']), flow.F0ElCp(p['M_inf']), phiInfFun))
        # add initial condition
        pbl.add(flow.Initial(self.msh, p['Fluid'], phiInfFun))
        # add farfield boundary conditions
        pbl.add(flow.Dirichlet(self.msh, p['Farfield'][0], phiInfFun))
        for i in range (1, len(p['Farfield'])):
            pbl.add(flow.Freestream(self.msh, p['Farfield'][i], velInfFun))
        # add solid boundaries and identify f/s boundary
        if p['Dim'] == 2:
            self.boundary = flow.Boundary(self.msh, [p['Wing'], p['Fluid']])
            pbl.add(self.boundary)
        else:
            for bd in p['Wings']:
                bnd = flow.Boundary(self.msh, [bd, p['Fluid']])
                pbl.add(bnd)
                if bd == p['Fsi']:
                    self.boundary = bnd
            if 'Fuselage' in p:
                pbl.add(flow.Boundary(self.msh, [p['Fuselage'], p['Fluid']]))
        # add wake/kutta condition
        if p['Dim'] == 2:
            pbl.add(flow.Wake(self.msh, [p['Wake'], p['Wake']+'_', p['Fluid']]))
            pbl.add(flow.Kutta(self.msh, [p['Te'], p['Wake']+'_', p['Wing'], p['Fluid']]))
        else:
            for i in range(0, len(p['Wakes'])):
                pbl.add(flow.Wake(self.msh, [p['Wakes'][i], p['Wakes'][i]+'_', p['Fluid'], p['TeTips'][i]]))

        # initialize the linear (inner) solver
        if p['LSolver'] == 'Pardiso':
            linsol = LinearSolver().pardiso()
        elif p['LSolver'] == 'GMRES':
            linsol = tbox.Gmres()
            linsol.setFillFactor(p['G_fill'])
            linsol.setRestart(p['G_restart'])
            linsol.setTolerance(p['G_tol'])
        elif p['LSolver'] == 'MUMPS':
            linsol = LinearSolver().mumps()
        elif p['LSolver'] == 'SparseLU':
            linsol = tbox.SparseLu()
        else:
            raise RuntimeError('Available linear solvers: Pardiso, GMRES, MUMPS or SparseLU, but ' + p['LSolver'] + ' was given!\n')
        # initialize the nonlinear (outer) solver
        if p['NSolver'] == 'Picard':
            self.solver = flow.Picard(pbl, linsol)
            self.solver.relax = p['Relaxation']
        elif p['NSolver'] == 'Newton':
            self.solver = flow.Newton(pbl, linsol)
            self.solver.lsTol = p['LS_tol']
            self.solver.maxLsIt = p['Max_it_LS']
            self.solver.avThrsh = p['AV_thrsh']
        else:
            raise RuntimeError('Available nonlinear solvers: Picard or Newton, but ' + p['NSolver'] + ' was given!\n')
        self.solver.nthreads = _nthreads
        self.solver.relTol = p['Rel_tol']
        self.solver.absTol = p['Abs_tol']
        self.solver.maxIt = p['Max_it']
        print('Linear solver: ', linsol)
        print("Number of threads: ", self.solver.nthreads)
        print("Maximum number of iterations: ", self.solver.maxIt)
        print("Objective relative residual: ", self.solver.relTol)
        print("Objective absolute residual: ", self.solver.absTol)
        print('\n')
        
    def run(self, t1, t2):
        """Run the solver for one steady (time) iteration.
        Adrien Crovato
        """
        #exeOK = self.solver.run()
        if not self.solver.run():
            raise RuntimeError('Flow solver diverged!\n')
        self.__setCurrentState()
    
    def __setCurrentState(self):
        """Compute nodal forces from nodal normalized forces
        Adrien Crovato
        """
        i = 0
        for n in self.boundary.nodes:
            self.nodalLoad_X[i] = -self.dynP * self.boundary.cLoadX[i]
            self.nodalLoad_Y[i] = -self.dynP * self.boundary.cLoadY[i]
            self.nodalLoad_Z[i] = -self.dynP * self.boundary.cLoadZ[i]
            i += 1

    def getNodalInitialPositions(self):
        """Get the initial position of each node
        Adrien Crovato
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
        Adrien Crovato
        """
        no = self.boundary.nodes[iVertex].no
        return no

    def applyNodalDisplacements(self, dx, dy, dz, dx_nM1, dy_nM1, dz_nM1, haloNodesDisplacements, time):
        """Apply displacements coming from solid solver to f/s interface after saving
        Adrien Crovato
        """
        self.mshDef.savePos()
        for i in range(self.boundary.nodes.size()):
            self.boundary.nodes[i].pos[0] = self.nodalInitPosX[i] + dx[i]
            self.boundary.nodes[i].pos[1] = self.nodalInitPosY[i] + dy[i]
            self.boundary.nodes[i].pos[2] = self.nodalInitPosZ[i] + dz[i]

    def meshUpdate(self, nt):
        """Deform the mesh using linear elasticity equations
        Adrien Crovato
        """
        self.mshDef.deform()
        
    def save(self, nt):
        """Save data on disk at each converged timestep
        Adrien Crovato
        """
        self.solver.save(nt, self.mshWriter)
        self.mshWriter.save(self.msh.name + "_" + str(nt))

    def initRealTimeData(self):
        """Initialize history file
        Adrien Crovato
        """
        histFile = open('FlowHistory.dat', 'w')
        histFile.write('{0:>12s}   {1:>12s}   {2:>12s}   {3:>12s}   {4:>12s}\n'.format('Time', 'FSI_Iter', 'C_Lift', 'C_Drag', 'C_Moment'))
        histFile.close()

    def saveRealTimeData(self, time, nFSIIter):
        """Save data at each fsi iteration
        Adrien Crovato
        """
        # history at each iteration
        histFile = open('FlowHistory.dat', 'a')
        histFile.write('{0:12.6f}   {1:12d}   {2:12.6f}   {3:12.6f}   {4:12.6f}\n'.format(time, nFSIIter, self.solver.Cl, self.solver.Cd, self.solver.Cm))
        histFile.close()
        # full solution at user-defined frequency
        if np.mod(nFSIIter+1, self.saveFreq) == 0:
            self.solver.save(1000000+old_div(int(nFSIIter+1),int(self.saveFreq)), self.mshWriter)

    def printRealTimeData(self, time, nFSIIter):
        """Print data on screen at the end of fsi simulation
        Adrien Crovato
        """
        print('[Flow lift, drag, moment]: {0:6.3f}, {1:6.4f}, {2:6.3f}'.format(self.solver.Cl, self.solver.Cd, self.solver.Cm))
        print('')
    
    def exit(self):
        """
        Exit the Flow solver
        """
        del self.solver
        del self.mshDef
        del self.msh
        del self.mshWriter
