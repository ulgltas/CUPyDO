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

@file FlowInterface.py
@brief Python interface between Flow and CUPyDO.
@author A. Crovato
"""

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import sys
import numpy as np
from cupydo.genericSolvers import FluidSolver

# ----------------------------------------------------------------------
#  FlowSolver class
# ----------------------------------------------------------------------

class Flow(FluidSolver):
    def __init__(self, _module, _nthreads, _saveFreq = sys.maxsize):
        # load the python module and initialize the solver
        module = __import__(_module)
        self.initFlow(module.getParams(), _nthreads)

        # count fsi nodes and get their positions
        self.nNodes = self.boundary.nodes.size()
        self.nHaloNode = 0
        self.nPhysicalNodes = self.nNodes - self.nHaloNode
        self.nodalInitPosX, self.nodalInitPosY, self.nodalInitPosZ = self.getNodalInitialPositions()

        # init save frequency (fsi)
        self.saveFreq = _saveFreq

        # generic init
        FluidSolver.__init__(self)

    def initFlow(self, p, _nthreads):
        """Initilize flow classes
        Adrien Crovato
        """
        import flow
        import tbox
        import tbox.gmsh as gmsh

        # basic checks
        if p['Dim'] != 2 and p['Dim'] != 3:
            raise Exception('Problem dimension should be 2 or 3, but ' + p['Dim'] + ' was given!\n')
        # basic config
        if p['Format'] == 'vtk':
            try:   
               import tboxVtk
               Writer = tboxVtk.VtkExport
               print "Found VTK libraries! Results will be saved in VTK format.\n"
            except:
               Writer = tbox.GmshExport
               print "VTK libraries not found! Results will be saved in gmsh format.\n"
        else:
            Writer = tbox.GmshExport
            print "Results will be saved in gmsh format.\n"

        # mesh the geometry
        self.msh = gmsh.MeshLoader(p['File'],__file__).execute(**p['Pars'])
        self.mshWriter = Writer(self.msh)
        gmshWriter = tbox.GmshExport(self.msh)
        if p['Dim'] == 2:
            mshCrck = tbox.MshCrack(self.msh, p['Dim'], gmshWriter, p['Wake'], [p['Fluid'], p['Fluid']+'_', p['Body'], p['Body']+'_', p['Farfield'][-1], p['Farfield'][-1]+'_'])
        else:
            mshCrck = tbox.MshCrack(self.msh, p['Dim'], gmshWriter, p['Wake'], [p['Fluid'], p['Fluid']+'_', p['Body'], p['Body']+'_', p['Farfield'][-1], p['Farfield'][-1]+'_', p['Symmetry'], p['Symmetry']+'_'], p['WakeTip'])
        del gmshWriter
        del mshCrck

        # initialize mesh deformation handler
        self.mshDef = tbox.MshDeform(self.msh, p['Dim'])
        self.mshDef.nthreads = _nthreads
        self.mshDef.setField(p['Fluid'])
        self.mshDef.setFixed(p['Farfield'])
        self.mshDef.setMoving([p['Body']])
        if p['Dim'] == 3:
            self.mshDef.setSymmetry([p['Symmetry']], 1)
        self.mshDef.setInternal([p['Wake'], p['Wake']+'_'])

        # initialize the problem
        p['AoA'] = p['AoA']*np.pi/180 # convert to radians
        phiInfFun = flow.Fun0PosPhiInf(p['Dim'], p['AoA'])
        if p['Dim'] == 2:
            velInfFun = tbox.Fct1C(-np.cos(p['AoA']), -np.sin(p['AoA']), 0.)
        else:
            velInfFun = tbox.Fct1C(-np.cos(p['AoA']), 0., -np.sin(p['AoA']))
        self.dynP = p['P_dyn']    
        pbl = flow.Problem(self.msh, p['Dim'], p['AoA'], p['M_inf'], p['S_ref'], p['c_ref'], p['x_ref'], p['y_ref'], p['z_ref'])

        # add medium
        if p['M_inf'] == 0:
            self.fCp = flow.Fun0EleCpL()
            pbl.set(flow.Medium(self.msh, p['Fluid'], flow.Fun0EleRhoL(), flow.Fun0EleMachL(), self.fCp, phiInfFun))
        else:
            self.fCp = flow.Fun0EleCp(p['gamma'], p['M_inf'])
            pbl.set(flow.Medium(self.msh, p['Fluid'], flow.Fun0EleRho(p['gamma'], p['M_inf'], p['M_crit']), flow.Fun0EleMach(p['gamma'], p['M_inf']), self.fCp, phiInfFun))
        # add initial condition
        pbl.add(flow.Assign(self.msh, p['Fluid'], phiInfFun), "IC")
        # add farfield and symmetry boundary conditions
        for bnd in p['Farfield']:
            pbl.add(flow.Neumann(self.msh, bnd, velInfFun))
        if p['Dim'] == 3:
            pbl.add(flow.Neumann(self.msh, p['Symmetry'], tbox.Fct1C(0., 0., 0.)))
        # add slip boundary condition and identify f/s boundary
        pbl.add(flow.Neumann(self.msh, p['Body'], tbox.Fct1C(0., 0., 0.)))
        self.boundary = flow.Boundary(self.msh, [p['Body'], p['Fluid']])
        pbl.add(self.boundary)
        # add wake/kutta condition
        if p['Dim'] == 2:
            pbl.add(flow.Wake(self.msh, [p['Wake'], p['Wake']+'_', p['Fluid']]))
            pbl.add(flow.Kutta(self.msh, [p['Te'], p['Wake']+'_', p['Body'], p['Fluid']]))
        else:
            pbl.add(flow.Wake(self.msh, [p['Wake'], p['Wake']+'_', p['Fluid'], p['TeTip']]))

        # initialize the solver
        if p['NSolver'] == 'Picard':
            self.solver = flow.Picard(pbl)
            self.solver.relax = p['Relaxation']
        elif p['NSolver'] == 'Newton':
            self.solver = flow.Newton(pbl)
            self.solver.lsTol = p['LS_tol']
            self.solver.maxLsIt = p['Max_it_LS']
            self.solver.avThrsh = p['AV_thrsh']
        else:
            raise RuntimeError('Available nonlinear solver type: Picard or Newton, but ' + p['NSolver'] + ' was given!\n')
        self.solver.nthreads = _nthreads
        self.solver.relTol = p['Rel_tol']
        self.solver.absTol = p['Abs_tol']
        self.solver.maxIt = p['Max_it']
        print "Number of threads: ", self.solver.nthreads
        print "Maximum number of iterations: ", self.solver.maxIt
        print "Objective relative residual: ", self.solver.relTol
        print "Objective absolute residual: ", self.solver.absTol
        print '\n'
        
    def run(self, t1, t2):
        """Run the solver for one steady (time) iteration.
        Adrien Crovato
        """
        #exeOK = self.solver.run()
        if not self.solver.run():
            raise RuntimeError('Flow solver diverged!\n')
        self.__setCurrentState()
    
    def __setCurrentState(self):
        """Compute nodal forces from nodal Pressure coefficient
        Adrien Crovato
        """
        # integrate Cp at element
        cpiE = self.boundary.integrate(self.solver.phi, self.fCp)
        # transfer integrated Cp from elements to nodes
        cfN = self.boundary.transfer(cpiE)
        i = 0
        for n in self.boundary.nodes:
            self.nodalLoad_X[i] = -self.dynP * cfN[i][0]
            self.nodalLoad_Y[i] = -self.dynP * cfN[i][1]
            self.nodalLoad_Z[i] = -self.dynP * cfN[i][2]
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
            x0[i] = n.pos.x[0]
            y0[i] = n.pos.x[1]
            z0[i] = n.pos.x[2]

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
            self.boundary.nodes[i].pos.x[0] = self.nodalInitPosX[i] + dx[i]
            self.boundary.nodes[i].pos.x[1] = self.nodalInitPosY[i] + dy[i]
            self.boundary.nodes[i].pos.x[2] = self.nodalInitPosZ[i] + dz[i]

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
        # History at each iteration
        histFile = open('FlowHistory.dat', 'a')
        histFile.write('{0:12.6f}   {1:12d}   {2:12.6f}   {3:12.6f}   {4:12.6f}\n'.format(time, nFSIIter, self.solver.Cl, self.solver.Cd, self.solver.Cm))
        histFile.close()
        # Full solution at user-defined frequency
        if np.mod(nFSIIter+1, self.saveFreq) == 0:
            self.solver.save(1000000+int(nFSIIter+1)/int(self.saveFreq), self.mshWriter)

    def printRealTimeData(self, time, nFSIIter):
        """Print data on screen at the end of fsi simulation
        Adrien Crovato
        """
        print '[Flow lift, drag, moment]: {0:6.3f}, {1:6.4f}, {2:6.3f}'.format(self.solver.Cl, self.solver.Cd, self.solver.Cm)
        print ''
    
    def exit(self):
        """
        Exit the Flow solver
        """
        del self.fCp
        del self.solver
        del self.mshDef
        del self.msh
        del self.mshWriter
