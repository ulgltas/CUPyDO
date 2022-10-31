from ..genericSolvers import FluidSolver
import pfem3Dw as w
import numpy as np

# %% Reads a Lua Input File

def readLua(path):

    with open(path,'r') as file: text = file.read()
    text = text.replace('"','').replace("'",'').replace(' ','')
    text = [x.split('=') for x in text.splitlines() if '=' in x]
    return dict(text)

# %% Interface Between PFEM3D and CUPyDO

class Pfem3D(FluidSolver):
    def __init__(self,param):

        path = param['cfdFile']
        input = readLua(path)

        # Read data from the lua file

        self.ID = input['Problem.id']
        self.group = input['Problem.interface']
        self.maxFactor = int(input['Problem.maxFactor'])
        self.autoRemesh = (input['Problem.autoRemeshing'] == 'true')

        # Problem class and functions initialization

        if self.ID == 'IncompNewtonNoT':

            self.run = self.runIncomp
            self.problem = w.ProbIncompNewton(path)
            self.applyDispBC = self.applyDispIncomp

        elif self.ID == 'WCompNewtonNoT':

            self.run = self.runWcomp
            self.problem = w.ProbWCompNewton(path)
            self.applyDispBC = self.applyDispWcomp

        else: raise Exception('Problem type not supported')

        # Stores the important objects and variables

        self.solver = self.problem.getSolver()
        self.prevSolution = w.SolutionData()
        self.mesh = self.problem.getMesh()
        self.dim = self.mesh.getDim()
        self.FSI = w.VectorInt()

        # FSI data and stores the previous time step 

        self.problem.copySolution(self.prevSolution)
        self.mesh.getNodesIndex(self.group,self.FSI)
        self.mesh.setComputeNormalCurvature(True)
        self.nPhysicalNodes = self.FSI.size()
        self.nNodes = self.nPhysicalNodes

        # Initializes the simulation data

        self.disp = np.zeros((self.nPhysicalNodes,3))
        self.initPos = self.getPosition()
        self.reload = False
        self.factor = 1
        self.ok = True

        # Prints the initial solution and stats

        self.problem.updateTime(-param['dt'])
        self.dt = param['dt']

        FluidSolver.__init__(self)
        self.problem.displayParams()
        self.problem.dump()

# %% Run for Incompressible Flows

    def runIncomp(self,t1,t2):

        print('\nSolve ({:.5e}, {:.5e})'.format(t1,t2))
        print('----------------------------------')

        # The line order is important here

        if not (self.reload and self.ok): self.factor //= 2
        self.factor = max(1,self.factor)
        self.resetSystem(t2-t1)
        iteration = 0

        # Main solving loop for the FSPC time step

        while iteration < self.factor:
            
            iteration += 1
            dt = (t2-t1)/self.factor
            self.solver.setTimeStep(dt)
            self.timeStats(dt+self.problem.getCurrentSimTime(),dt)
            self.ok = self.solver.solveOneTimeStep()

            if not self.ok:

                print('PFEM3D: Problem occured\n')
                if 2*self.factor > self.maxFactor: return False
                self.factor = 2*self.factor
                self.resetSystem(t2-t1)
                iteration = 0

        print('PFEM3D: Successful run')
        self.setCurrentState()
        return True
        
# %% Run for Weakly Compressible Flows

    def runWcomp(self,t1,t2):

        print('\nSolve ({:.5e}, {:.5e})'.format(t1,t2))
        print('----------------------------------')

        # Estimate the time step only once
        
        self.resetSystem(t2-t1)
        self.solver.computeNextDT()
        self.factor = int((t2-t1)/self.solver.getTimeStep())
        if self.factor > self.maxFactor: return False
        dt = (t2-t1)/self.factor
        self.timeStats(t2,dt)
        iteration = 0

        # Main solving loop for the FSPC time step

        while iteration < self.factor:
    
            iteration += 1
            self.solver.setTimeStep(dt)
            self.solver.solveOneTimeStep()

        print('PFEM3D: Successful run')
        self.setCurrentState()
        return True

# %% Apply Boundary Conditions

    def applyNodalDisplacements(self,dx,dy,dz,*_):
        self.disp = np.transpose([dx,dy,dz])

    # For implicit and incompressible flows

    def applyDispIncomp(self,distance,dt):

        BC = (distance)/dt
        for i in range(self.dim):
            for j,k in enumerate(self.FSI):
                self.mesh.setNodeState(k,i,BC[j,i])

    # For explicit weakly compressive flows

    def applyDispWcomp(self,distance,dt):

        velocity = self.getVelocity()
        BC = 2*(distance-velocity*dt)/(dt*dt)

        # Update the FSI node states BC

        for i in range(self.dim):
            idx = int(self.dim+2+i)

            for j,k in enumerate(self.FSI):
                self.mesh.setNodeState(k,idx,BC[j,i])

# %% Return Nodal Values

    def getPosition(self):

        pos = self.disp.copy()
        for i in range(self.dim):
            for j,k in enumerate(self.FSI):
                pos[j,i] = self.mesh.getNode(k).getCoordinate(i)

        return pos

    # Computes the nodal velocity vector

    def getVelocity(self):

        vel = self.disp.copy()
        for i in range(self.dim):
            for j,k in enumerate(self.FSI):
                vel[j,i] = self.mesh.getNode(k).getState(i)

        return vel

    # Computes the reaction nodal loads

    def setCurrentState(self):

        vec = w.VectorArrayDouble3()
        self.solver.computeLoads(self.group,self.FSI,vec)

        for i in range(self.nNodes):

            self.nodalLoad_X[i] = -vec[i][0]
            self.nodalLoad_Y[i] = -vec[i][1]
            self.nodalLoad_Z[i] = -vec[i][2]

# %% Other Functions

    def update(self,dt):

        self.mesh.remesh(False)
        if self.ID == 'IncompNewtonNoT': self.solver.precomputeMatrix()
        self.problem.copySolution(self.prevSolution)
        FluidSolver.update(self,dt)
        self.reload = False

    # Prepare to solve one time step

    def resetSystem(self,dt):

        if self.reload: self.problem.loadSolution(self.prevSolution)
        if self.autoRemesh and (self.ID == 'IncompNewtonNoT'):
            if self.reload: self.solver.precomputeMatrix()

        distance = self.disp-(self.getPosition()-self.initPos)
        self.applyDispBC(distance,dt)
        self.reload = True

    # Other utilitary functions

    def getNodalIndex(self,index):
        return index

    def getNodalInitialPositions(self):
        return np.transpose(self.initPos)

# %% Print Functions

    def exit(self):

        print('======================================')
        self.problem.displayTimeStats()
        print('======================================')
        print('\nExit PFEM3D')

    # Save te results into a file

    def save(self,_):
        self.problem.dump()

    # Display the current simulation state

    def timeStats(self,time,dt):

        start = self.problem.getCurrentSimTime()
        print('t1 = {:.5e} - dt = {:.3e}'.format(start,dt))
        print('t2 = {:.5e} - factor = {:.0f}'.format(time,self.factor))
        print('----------------------------------')
