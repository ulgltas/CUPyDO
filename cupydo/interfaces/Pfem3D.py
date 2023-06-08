from ..genericSolvers import FluidSolver
from ..utilities import titlePrint
import pfem3Dw as w
import numpy as np

# %% Interface Between PFEM3D and CUPyDO

class Pfem3D(FluidSolver):
    def __init__(self,param):

        titlePrint('Initializing PFEM3D')
        self.problem = w.getProblem(param['cfdFile'])

        # Incompressible or weakly compressible solver

        if 'WC' in self.problem.getID():
            
            self.implicit = False
            self.run = self.runExplicit
            self.maxDivision = 2000

        else:
            
            self.implicit = True
            self.run = self.runImplicit
            self.maxDivision = 10

        # Stores the important objects and variables

        self.FSI = w.VectorInt()
        self.mesh = self.problem.getMesh()
        self.mesh.getNodesIndex('FSInterface',self.FSI)
        self.solver = self.problem.getSolver()
        self.nPhysicalNodes = self.FSI.size()
        self.nNodes = self.FSI.size()

        # Initialize the boundary conditions

        self.BC = list()
        self.dim = self.mesh.getDim()

        for i in self.FSI:

            vector = w.VectorDouble(3)
            self.mesh.getNode(i).setExtState(vector)
            self.BC.append(vector)

        # Save mesh after initializing the BC pointer

        self.prevSolution = w.SolutionData()
        self.problem.copySolution(self.prevSolution)
        self.problem.displayParams()
        self.problem.dump()

        # Store temporary simulation variables

        self.disp = np.zeros((self.nPhysicalNodes,3))
        self.initPos = self.getPosition()
        self.vel = self.getVelocity()
        self.dt = param['dt']
        
        FluidSolver.__init__(self)

# %% Run for implicit integration scheme

    def runImplicit(self,t1,t2):

        print('\nt = {:.5e} - dt = {:.5e}'.format(t2,t2-t1))
        self.problem.loadSolution(self.prevSolution)
        dt = float(t2-t1)
        count = int(1)

        # Main solving loop for the fluid simulation

        while count > 0:
            
            self.solver.setTimeStep(dt)
            if not self.solver.solveOneTimeStep():
                
                dt = float(dt/2)
                count = np.multiply(2,count)
                if dt < (t2-t1)/self.maxDivision:
                    raise Exception('Too large time step')
                continue

            count = count-1
        self.__setCurrentState()
        return True

# %% Run for explicit integration scheme

    def runExplicit(self,t1,t2):

        print('\nt = {:.5e} - dt = {:.5e}'.format(t2,t2-t1))
        self.problem.loadSolution(self.prevSolution)
        iteration = 0

        # Estimate the time step for stability

        self.solver.computeNextDT()
        division = int((t2-t1)/self.solver.getTimeStep())
        if division > self.maxDivision:
            raise Exception('Too large time step')
        dt = (t2-t1)/division

        # Main solving loop for the fluid simulation

        while iteration < division:
    
            iteration += 1
            self.solver.setTimeStep(dt)
            self.solver.solveOneTimeStep()

        self.__setCurrentState()
        return True

# %% Apply Boundary Conditions

    def applyNodalDisplacements(self,dx,dy,dz,*_):

        BC = (np.transpose([dx,dy,dz])-self.disp)/self.dt
        if not self.implicit: BC = 2*(BC-self.vel)/self.dt

        for i,vector in enumerate(BC):
            for j,val in enumerate(vector): self.BC[i][j] = val

# %% Return Nodal Values

    def getPosition(self):

        vector = np.zeros(self.disp.shape)

        for i in range(self.dim):
            for j,k in enumerate(self.FSI):
                vector[j,i] = self.mesh.getNode(k).getCoordinate(i)

        return vector

    # Computes the nodal velocity vector

    def getVelocity(self):

        vector = np.zeros(self.disp.shape)
        
        for i in range(self.dim):
            for j,k in enumerate(self.FSI):
                vector[j,i] = self.mesh.getNode(k).getState(i)

        return vector

    # Computes the reaction nodal loads

    def __setCurrentState(self):

        vector = w.VectorVectorDouble()
        self.solver.computeLoads('FSInterface',self.FSI,vector)

        for i in range(self.nNodes):

            self.nodalLoad_X[i] = -vector[i][0]
            self.nodalLoad_Y[i] = -vector[i][1]
            if self.dim == 3: self.nodalLoad_Z[i] = -vector[i][2]

# %% Other Functions

    def update(self,_):

        self.mesh.remesh(False)
        if self.implicit: self.solver.precomputeMatrix()
        self.problem.copySolution(self.prevSolution)
        self.disp = self.getPosition()-self.initPos
        self.vel = self.getVelocity()

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
        print('t2 = {:.5e} - division = {:.0f}'.format(time,self.factor))
        print('----------------------------------')
