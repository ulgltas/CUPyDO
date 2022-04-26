from ..genericSolvers import FluidSolver
from slpp import slpp as lua
from sys import stdout
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

        print('\nInitializing PFEM3D')
        path = param['cfdFile']
        input = read(path)

        # Read some input parameters

        self.ID = input['id']
        self.group = input['FSInterface']

        # Problem class initialization

        if self.ID == 'IncompNewtonNoT': self.problem = w.ProbIncompNewton(path)
        elif self.ID == 'WCompNewtonNoT': self.problem = w.ProbWCompNewton(path)
        else: raise Exception('Problem type not supported')

        # Gets some important objects and variables

        self.solver = self.problem.getSolver()
        self.prevSolution = w.SolutionData()
        self.mesh = self.problem.getMesh()
        self.load = w.VectorArrayDouble3()
        self.FSI = w.VectorInt()

        # Set the initial time step

        self.problem.updateTime(-param['dt'])
        self.solver.setTimeStep(1e16)
        self.solver.computeNextDT()
        self.dt = param['dt']

        # Number of nodes at the FSInterface
    
        self.mesh.getNodesIndex(self.group,self.FSI)
        self.nPhysicalNodes = self.FSI.size()
        self.nNodes = self.FSI.size()
        self.dim = self.mesh.getDim()

        # Initializes the simulation data

        self.pos0 = self.getPosition()
        self.problem.copySolution(self.prevSolution)
        self.prevDisp = np.zeros((self.nNodes,3))
        self.reload = False

        # Prints and initializes the solver

        FluidSolver.__init__(self)
        self.problem.displayParams()
        
# %% Computes one time increment

    def run(self,t1,t2):
        self.reload = True

        while True:

            dt = self.solver.getTimeStep()
            time = self.problem.getCurrentSimTime()

            if (dt+time-t2)/dt > -0.2:
                
                self.solver.setTimeStep(t2-time)
                ok = self.solver.solveOneTimeStep()
                self.printSolverTime(ok)
                self.solver.computeNextDT()
                if ok: break

            else:
                ok = self.solver.solveOneTimeStep()
                self.printSolverTime(ok)
                self.solver.computeNextDT()

        # Computes nodal fluid loads

        self.setCurrentState()
        return True

    # Prints and set next time step

    def printSolverTime(self,ok):

        dt = self.solver.getTimeStep()
        time = self.problem.getCurrentSimTime()

        if self.ID == 'IncompNewtonNoT':
            print('PFEM3D [{:<1}] : t = {:.5e} - dt = {:.3e}'.format(ok,time,dt))

        elif self.ID == 'WCompNewtonNoT':
            stdout.write('\rPFEM3D : t = {:.5e} - dt = {:.3e}'.format(time,dt))
            stdout.flush()

# %% Get and Set Nodal Values

    def getPosition(self):
        pos = np.zeros((self.nNodes,3))

        for i in range(self.nNodes):
            node = self.mesh.getNode(self.FSI[i])
            pos[i] = [node.getCoordinate(j) for j in range(3)]

        return pos

    def getVelocity(self):
        vel = np.zeros((self.nNodes,3))

        for i in range(self.nNodes):
            node = self.mesh.getNode(self.FSI[i])
            for j in range(self.dim): vel[i,j] = node.getState(j)

        return vel

    # Sets current nodal loads on the FSI

    def setCurrentState(self):

        self.solver.computeLoads(self.group,self.FSI,self.load)
        for i in range(self.nNodes):

            self.nodalLoad_X[i] = -self.load[i][0]
            self.nodalLoad_Y[i] = -self.load[i][1]
            self.nodalLoad_Z[i] = -self.load[i][2]

    # Some other utilitary functions

    def getNodalIndex(self,index):
        return self.FSI[index]

    def getNodalInitialPositions(self):
        return np.transpose(self.pos0)

# %% Sets Boundary Conditions

    def applyNodalDisplacements(self,dx,dy,dz,dx_nM1,dy_nM1,dz_nM1,haloNodesDisplacements,time):
        
        dDisp = np.transpose([dx,dy,dz])-self.prevDisp
        if self.reload: self.problem.loadSolution(self.prevSolution)

        # Computes the state according to Metafor

        if (self.ID == 'IncompNewtonNoT'):

            index = int(0)
            BC = dDisp/self.dt

        elif (self.ID == 'WCompNewtonNoT'):

            index = int(self.dim+2)
            vel = self.getVelocity()
            BC = 2*(dDisp-vel*self.dt)/(self.dt*self.dt)

        # Update the FSI node states BC

        for i in range(self.nNodes):
            for j in range(self.dim):
                self.mesh.setNodeState(self.FSI[i],index+j,BC[i,j])

# %% Update and Save Results

    def update(self,dt):

        self.mesh.remesh(False)
        if (self.ID == 'IncompNewtonNoT'): self.solver.precomputeMatrix()
        self.prevDisp = self.getPosition()-self.pos0
        self.problem.copySolution(self.prevSolution)
        FluidSolver.update(self,dt)
        self.reload = False

    def save(self,nt):
        self.problem.dump()

# %% Exits the Solver

    def exit(self):

        print('======================================')
        self.problem.displayTimeStats()
        print('======================================')
        print('\nExit PFEM3D')