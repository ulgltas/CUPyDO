from ..genericSolvers import FluidSolver
from ..utilities import titlePrint
import pfem3Dw as w
import numpy as np

# %% Reads a Lua Input File

def read(path):

    with open(path,'r') as file: text = file.read()
    text = text.replace('"','').replace("'",'').replace(' ','')
    text = [x.split('=') for x in text.splitlines() if '=' in x]
    return dict(text)

# %% Interface Between PFEM3D and CUPyDO

class Pfem3D(FluidSolver):
    def __init__(self,param):

        titlePrint('Initializing PFEM3D')

        path = param['cfdFile']
        self.problem = w.getProblem(path)
        self.autoRemesh = self.problem.hasAutoRemeshing()
        problemID = self.problem.getID()

        # Read the input Lua file

        input = read(path)
        self.group = input['Problem.interface']
        self.maxFactor = int(input['Problem.maxFactor'])
        if problemID[:2] == 'WC': self.implicit = False
        else: self.implicit = True

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
        self.problem.displayParams()

        # Initializes the simulation data

        self.disp = np.zeros((self.nPhysicalNodes,3))
        self.initPos =  self.getPosition()
        self.dt = param['dt']
        self.reload = False
        self.factor = 1
        self.ok = True
        
        FluidSolver.__init__(self)

# %% Calculates One Time Step

    def run(self,t1,t2):

        print('\nSolve ({:.5e}, {:.5e})'.format(t1,t2))
        print('----------------------------------')
        if self.implicit: return self.runImplicit(t1,t2)
        else: return self.runExplicit(t1,t2)

    # Run for implicit integration scheme

    def runImplicit(self,t1,t2):

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

        self.__setCurrentState()
        return True

    # Run for explicit integration scheme

    def runExplicit(self,t1,t2):
        
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

        self.__setCurrentState()
        return True

# %% Apply Boundary Conditions

    def applyNodalDisplacements(self,dx,dy,dz,*_):
        self.disp = np.transpose([dx,dy,dz])

    # Update and apply the nodal displacement

    def applyDispBC(self,distance,dt):

        if self.implicit:

            BC = w.VectorVectorDouble(distance/dt)
            self.solver.setVelocity(self.FSI,BC)

        else:

            BC = 2*(distance-self.getVelocity()*dt)
            BC = w.VectorVectorDouble(BC/np.square(dt))
            self.solver.setAcceleration(self.FSI,BC)

# %% Return Nodal Values

    def getPosition(self):

        vector = np.zeros((self.nPhysicalNodes,3))

        for i in range(self.dim):
            for j,k in enumerate(self.FSI):
                vector[j,i] = self.mesh.getNode(k).getCoordinate(i)

        return vector

    # Computes the nodal velocity vector

    def getVelocity(self):

        vector = np.zeros((self.nPhysicalNodes,3))
        
        for i in range(self.dim):
            for j,k in enumerate(self.FSI):
                vector[j,i] = self.mesh.getNode(k).getState(i)

        return vector
        
    # Computes the reaction nodal loads

    def __setCurrentState(self):

        vector = w.VectorVectorDouble()
        self.solver.computeLoads(self.group,self.FSI,vector)

        for i in range(self.nNodes):

            self.nodalLoad_X[i] = -vector[i][0]
            self.nodalLoad_Y[i] = -vector[i][1]
            if self.dim == 3: self.nodalLoad_Z[i] = -vector[i][2]

# %% Other Functions

    def update(self,_):

        self.mesh.remesh(False)
        if self.implicit: self.solver.precomputeMatrix()
        self.problem.copySolution(self.prevSolution)
        self.reload = False

    # Prepare to solve one time step

    def resetSystem(self,dt):

        if self.reload: self.problem.loadSolution(self.prevSolution)
        if self.autoRemesh and self.implicit and self.reload:
            self.solver.precomputeMatrix()

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
