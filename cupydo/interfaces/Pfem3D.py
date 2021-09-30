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

        # Initializes the loads and displacements

        self.load = wraper.VectorArrayDouble3(1)
        self.dx = np.zeros(self.nPhysicalNodes)
        self.dy = np.zeros(self.nPhysicalNodes)
        self.dz = np.zeros(self.nPhysicalNodes)

        # Checks if autoRemeshing is false

        if self.solver.getAutoRemeshing():
            raise Exception('Auto remeshing must be disabled')

        # Initializes the fluid solver

        self.saveNode = self.mesh.copyNodesList()
        FluidSolver.__init__(self)
        self.initRealTimeData()
        
# %% Calculates One Increment From t1 to t2

    def run(self,t1,t2):

        self.mesh.setNodesList(self.saveNode)
        self.applyNodalDisplacements(self.dx,self.dy,self.dz,0,0,0,0,0)

        # Initialization of the simulation
        
        self.problem.updateTime(t1-t2)
        time = self.problem.getCurrentSimTime()
        dt = self.solver.getTimeStep()
        progress = True
        print()

        # Solves without remeshing until t2

        while progress:
            if (dt+time-t2)/dt > -0.1:
                    
                self.solver.forceTimeStep(t2-time)
                ok = self.solver.solveOneTimeStep()
                time = self.problem.getCurrentSimTime()
                dt = self.solver.getTimeStep()
                self.solver.computeNextDT()
                if(ok): progress = False

            else:
                ok = self.solver.solveOneTimeStep()
                time = self.problem.getCurrentSimTime()
                dt = self.solver.getTimeStep()
                self.solver.computeNextDT()

            # Prints the current state
            
            print(ok,': t = {:.6f} - dt = {:.2e}'.format(time,dt))
            dt = self.solver.getTimeStep()

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

        self.solver.computeLoads(self.load)
        self.nodalLoad_X = np.zeros(self.nPhysicalNodes)
        self.nodalLoad_Y = np.zeros(self.nPhysicalNodes)
        self.nodalLoad_Z = np.zeros(self.nPhysicalNodes)
        
        for i in range(self.nPhysicalNodes):
            
            idx = self.nodeID[i]
            self.nodalLoad_X[i] = self.load[idx][0]
            self.nodalLoad_Y[i] = self.load[idx][1]
            self.nodalLoad_Z[i] = self.load[idx][2]

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
        self.saveNode = self.mesh.copyNodesList()
        self.updateID()

        self.problem.updateTime(dt)
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

        self.problem.writeExtractors()
        return

    def initRealTimeData(self):

        self.problem.dump()
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