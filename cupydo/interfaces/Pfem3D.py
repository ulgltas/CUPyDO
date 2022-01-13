from ..genericSolvers import FluidSolver
from slpp import slpp as lua
import pfem3Dw as w
import numpy as np
import os

# %% Translate the Lua Table Into a Dict

def read(path):

    if(os.path.isfile(path)):
        with open(path,'r') as file:

            input = file.read().replace(' ','')
            input = input.replace('Problem=','')
            input = lua.decode(input)

    else: raise Exception('Cannot open lua file '+path)
    return input

# %% Interface Between PFEM3D and CUPyDO

class Pfem3D(FluidSolver):
    def __init__(self,param):

        print('\n***************************** Initializing PFEM3D *****************************')
        path = param['cfdFile']
        input = read(path)

        # Problem class initialization

        if input['id'] == 'IncompNewtonNoT':

            self.problem = w.ProbIncompNewton(path)
            self.typeBC = 'velocity'

        elif input['id'] == 'WCompNewtonNoT':

            self.problem = w.ProbWCompNewton(path)
            self.typeBC = 'acceleration'

        else: raise Exception('Problem type not supported')
        if not input['useCupydo']: raise Exception('UseCupydo must be True')
        self.group = input['FSInterface']

        # Gets some important objects and variables

        self.solver = self.problem.getSolver()
        self.mesh = self.problem.getMesh()
        self.load = w.VectorArrayDouble3()
        self.prevMesh = w.Mesh()
        self.FSI = w.VectorInt()

        # Number of nodes at the FSInterface
    
        self.mesh.getNodesIndexOfTag(self.group,self.FSI)
        self.nPhysicalNodes = self.FSI.size()
        self.nNodes = self.FSI.size()
        self.dim = self.mesh.getDim()

        # Initializes the simulation data

        self.prevTime = self.problem.getCurrentSimTime()
        self.pos = self.getNodalInitialPositions()
        self.disp = np.zeros((self.nNodes,3))
        self.vel = np.zeros((self.nNodes,3))
        self.BC = np.zeros((self.nNodes,3))
        self.prevMesh.deepCopy(self.mesh)
        self.reload = False
        self.setNodalBC()

        # Prints and initializes the solver

        FluidSolver.__init__(self)
        self.problem.updateTime(-param['dt'])
        self.solver.setTimeStep(param['dt'])
        self.problem.displayParams()
        self.dt = param['dt']
        
# %% Calculates One Increment From t1 to t2

    def run(self,t1,t2):
        print('PFEM3D: Run from {:.5e} to {:.5e}'.format(t1,t2))

        # Solves with remeshing until t2

        while True:
            
            self.solver.computeNextDT()
            dt = self.solver.getTimeStep()
            time = self.problem.getCurrentSimTime()

            if (dt+time-t2)/dt > -0.2:
                
                dt = t2-time
                self.solver.setTimeStep(dt)
                if self.solver.solveOneTimeStep(): self.remesh(); break

            else:
                if self.solver.solveOneTimeStep(): self.remesh()

        # Computes nodal fluid loads

        self.setCurrentState()
        self.reload = True

# %% Get and Set Nodal Values

    def getNodalInitialPositions(self):

        coord = np.zeros((self.nNodes,3))
        for i in range(self.nNodes):

            idx = self.FSI[i]
            node = self.mesh.getNode(idx)
            coord[i] = [node.getCoordinate(j) for j in range(3)]

        return np.transpose(coord)

    # Sets current node states

    def setCurrentState(self):

        self.solver.computeLoads(self.group,self.FSI,self.load)
        for i in range(self.nNodes):

            self.nodalLoad_X[i] = -self.load[i][0]
            self.nodalLoad_Y[i] = -self.load[i][1]
            self.nodalLoad_Z[i] = -self.load[i][2]

    # Returns the index of the index-th interface node

    def getNodalIndex(self,index):
        return self.FSI[index]

    # Prescribes nodal positions from solid solver

    def applyNodalDisplacements(self,dx,dy,dz,dx_nM1,dy_nM1,dz_nM1,haloNodesDisplacements,time):

        if self.reload:
            self.problem.loadMesh(self.prevMesh,self.prevTime)
            self.mesh.getNodesIndexOfTag(self.group,self.FSI)

        # BC for incompressible fluids

        if self.typeBC == 'velocity':

            disp = np.transpose([dx,dy,dz])-self.disp
            self.BC = disp/self.dt
            self.setNodalBC()

        # BC for weakly compressible fluids

        elif self.typeBC == 'acceleration':

            for i in range(self.nNodes):
                for j in range(self.dim):
                    self.vel[i,j] = self.mesh.getNode(self.FSI[i]).getState(j)

            disp = np.transpose([dx,dy,dz])-self.disp
            self.BC = 2*(disp-self.vel*self.dt)/(self.dt*self.dt)
            self.setNodalBC()

# %% Update and Save Results

    def update(self,dt):

        self.prevMesh.deepCopy(self.mesh)
        self.prevTime = self.problem.getCurrentSimTime()
        self.disp = np.transpose(self.getNodalInitialPositions()-self.pos)
        FluidSolver.update(self,dt)
        self.reload = False

    def save(self,nt):
        self.problem.dump()

    def setNodalBC(self):

        nodalBC = w.MapIntArrayDouble3()
        for i in range(self.nNodes): nodalBC[self.FSI[i]] = self.BC[i]
        self.solver.setNodalBC(nodalBC)

    def remesh(self):
        
        self.mesh.remesh(False)
        self.mesh.getNodesIndexOfTag(self.group,self.FSI)
        self.setNodalBC()

# %% Exits the Solver

    def exit(self):

        print('======================================')
        self.problem.displayTimeStats()
        print('======================================')
        print('\n***************************** Exit PFEM3D *****************************')