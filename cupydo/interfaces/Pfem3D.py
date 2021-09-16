from ..genericSolvers import FluidSolver
from slpp import slpp as lua
import pfem3Dw as wraper
import numpy as np
import os

# %% Interface Between PFEM3D and CUPyDO

class Pfem3D(FluidSolver):
    def __init__(self,param):

        name = param['cfdFile']+'.lua'
        print(f"\nInitializing PFEM3D\n")

        # Reads the Lua table into a dict

        if(os.path.isfile(name)):
            with open(name,'r') as file:

                input = file.read().replace(' ','')
                input = input.replace('Problem=','')
                input = lua.decode(input)

        else: raise Exception('Cannot open lua file '+name)

        # Problem type availables

        incomp = ['IncompNewtonNoT','Bingham','Boussinesq','Conduction']
        wecomp = ['WCompNewtonNoT','BoussinesqWC']

        # Problem initialization

        if input['id'] in incomp: self.problem = wraper.ProbIncompNewton(name)
        elif input['id'] in wecomp: self.problem = wraper.ProbWCompNewton(name)
        else: raise Exception('Unknown problem type '+input['id'])

        # Gets some important objects and variables

        self.solver = self.problem.getSolver()
        self.mesh = self.problem.getMesh()

        # Ghost and physical nodes at the FS interface

        self.updateID()
        self.nNodes = len(self.nodeID)
        self.nPhysicalNodes = self.nNodes

        # Initializes the fluid solver

        self.dx = np.zeros(self.nPhysicalNodes)
        self.dy = np.zeros(self.nPhysicalNodes)
        self.dz = np.zeros(self.nPhysicalNodes)
        self.solver.forceTimeStep(param['dt'])
        FluidSolver.__init__(self)
        self.initRealTimeData()

# %% Calculates One Increment From t1 to t2

    def run(self,t1,t2):

        if(np.allclose(t1,self.time)):
            self.mesh.restoreNodesList()
            self.applyNodalDisplacements(self.dx,self.dy,self.dz,0,0,0,0,0)

        # Solve and reverse the simulation time
        
        self.time = t1
        self.problem.updateTime(t1-t2)
        ok = self.solver.solveOneTimeStep()
        if(not ok): raise Exception('Failed to solve time step')

        # Computes nodal fluid loads

        self.setCurrentState()
        return

# %% Get and Set Nodal Values

    def getNodalInitialPositions(self):

        coord = np.zeros((3,self.nPhysicalNodes))
        for i in range(self.nPhysicalNodes):

            idx = self.nodeID[i]
            node = self.mesh.getNode(idx)
            coord[0,i] = node.getCoordinate(0)
            coord[1,i] = node.getCoordinate(1)
            coord[2,i] = node.getCoordinate(2)

        return coord[0],coord[1],coord[2]

    # Sets current node states

    def setCurrentState(self):

        load = self.solver.computeLoads()
        self.nodalLoad_X = np.zeros(self.nPhysicalNodes)
        self.nodalLoad_Y = np.zeros(self.nPhysicalNodes)
        self.nodalLoad_Z = np.zeros(self.nPhysicalNodes)
        
        for i in range(self.nPhysicalNodes):
            
            idx = self.nodeID[i]
            self.nodalLoad_X[i] = load[idx][0]
            self.nodalLoad_Y[i] = load[idx][1]
            self.nodalLoad_Z[i] = load[idx][2]

    # Returns the index of the iVertex-th interface node

    def getNodalIndex(self,iVertex):
        return self.nodeID[iVertex]

    # Prescribes Nodal Positions from Solid Solver

    def applyNodalDisplacements(self,dx,dy,dz,dx_nM1,dy_nM1,dz_nM1,haloNodesDisplacements,time):

        for i in range(self.nPhysicalNodes):

            idx = self.nodeID[i]
            self.mesh.setNodePosition(self.dx0[i]+dx[i],idx,0)
            self.mesh.setNodePosition(self.dy0[i]+dy[i],idx,1)
            self.mesh.setNodePosition(self.dz0[i]+dz[i],idx,2)

        self.dx = dx
        self.dy = dy
        self.dz = dz
        return

# %% Update and Save Results

    def update(self,dt):

        self.mesh.remesh(self.problem.isOutputVerbose())
        self.mesh.saveNodesList()
        self.updateID()

        self.problem.updateTime(dt)
        self.solver.forceTimeStep(dt)
        FluidSolver.update(self,dt)
        return

    # Updates the FSI node index

    def updateID(self):
        
        self.nodeID = []
        for i in range(self.mesh.getNodesCount()):
            if self.mesh.getNodeType(i)=='FSInterface':
                self.nodeID.append(i)

# %% Prints and Save Temporary Outputs

    def save(self,nt):

        self.problem.dump(False)
        return

    def initRealTimeData(self):

        self.problem.dump(False)
        self.problem.displayParams()
        dx,dy,dz = self.getNodalInitialPositions()

        self.dx0 = dx.copy()
        self.dy0 = dy.copy()
        self.dz0 = dz.copy()
        self.time = -np.inf
        return

# %% Exits the Solver

    def exit(self):

        print('======================================')
        self.problem.displayTimeStats()
        print('======================================')
        print('\nExit PFEM3D')
        return