from ..genericSolvers import FluidSolver
import pfem3Dw as w
import numpy as np

# %% Interface Between PFEM3D and CUPyDO

class Pfem3D(FluidSolver):
    def __init__(self,param):

        print('\nInitializing PFEM3D')
        path = param['cfdFile']
        self.read(path)

        # Problem class and functions initialization

        if self.ID == 'IncompNewtonNoT':

            self.run = self.runIncomp
            self.problem = w.ProbIncompNewton(path)
            self.applyNodalDisplacements = self.applyDispIncomp

        elif self.ID == 'WCompNewtonNoT':

            self.run = self.runWcomp
            self.problem = w.ProbWCompNewton(path)
            self.applyNodalDisplacements = self.applyDispWcomp

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
        self.nNodes = self.FSI.size()

        # Set the initial time step

        self.problem.updateTime(-param['dt'])
        self.dt = param['dt']

        # Initializes the simulation data

        self.prevDisp = np.zeros((self.nNodes,3))
        self.disp = np.zeros((self.nNodes,3))
        self.load = np.zeros((self.nNodes,3))
        self.vel = np.zeros((self.nNodes,3))
        self.pos = np.zeros((self.nNodes,3))

        # Initializes some FSPC data

        self.reload = False
        self.factor = 1
        self.ok = True
        self.div = 2

        # Prints the initial solution and stats

        FluidSolver.__init__(self)
        self.problem.displayParams()
        self.problem.dump()

# %% Computes a Time Increment

    def runIncomp(self,t1,t2):

        if not (self.reload and self.ok): self.factor //= self.div
        self.factor = max(1,self.factor)
        self.reload = True
        iter = 0

        # Main solving loop for the CUPyDO time step

        while iter < self.factor:

            iter += 1
            dt = (t2-t1)/self.factor
            self.solver.setTimeStep(dt)
            self.timeStats(dt+self.problem.getCurrentSimTime(),dt)
            self.ok = self.solver.solveOneTimeStep()

            if not self.ok:

                self.applyNodalDisplacements(*self.disp.T)
                self.factor *= self.div
                iter = 0

        self.setCurrentState()
        return True

    # For explicit weakly compressive flows

    def runWcomp(self,t1,t2):

        iter = 0
        self.reload = True
        self.solver.computeNextDT()
        self.factor = int((t2-t1)/self.solver.getTimeStep())
        dt = (t2-t1)/self.factor
        self.timeStats(t2,dt)

        # Main solving loop for the CUPyDO time step

        while iter < self.factor:
    
            iter += 1
            self.solver.setTimeStep(dt)
            self.solver.solveOneTimeStep()

        self.setCurrentState()
        return True

# %% Sets Boundary Conditions

    def applyDispIncomp(self,dx,dy,dz,*_):

        disp = np.transpose([dx,dy,dz])
        if self.reload: self.problem.loadSolution(self.prevSolution)
        BC = (disp-self.prevDisp)/self.dt
        self.disp = disp.copy()

        # Update the FSI node states BC

        for i in range(self.nNodes):
            for j in range(self.dim):
                self.mesh.setNodeState(self.FSI[i],j,BC[i,j])

    # For explicit weakly compressive flows

    def applyDispWcomp(self,dx,dy,dz,*_):

        disp = np.transpose([dx,dy,dz])
        if self.reload: self.problem.loadSolution(self.prevSolution)
        BC = 2*(disp-self.prevDisp-self.getVelocity()*self.dt)/self.dt**2
        self.disp = disp.copy()
        idx = int(self.dim+2)

        # Update the FSI node states BC

        for i in range(self.nNodes):
            for j in range(self.dim):
                self.mesh.setNodeState(self.FSI[i],idx+j,BC[i,j])

# %% Gets Nodal Values

    def getPosition(self):

        for i in range(self.nNodes):
            node = self.mesh.getNode(self.FSI[i])
            self.pos[i] = [node.getCoordinate(j) for j in range(3)]

        return self.pos.copy()

    def getVelocity(self):

        for i in range(self.nNodes):
            node = self.mesh.getNode(self.FSI[i])
            self.vel[i] = [node.getState(j) for j in range(3)]

        return self.vel.copy()

    # Sets current nodal loads on the FSI

    def setCurrentState(self):

        load = w.VectorArrayDouble3()
        self.solver.computeLoads(self.group,self.FSI,load)
        for i in range(self.nNodes):

            self.nodalLoad_X[i] = -load[i][0]
            self.nodalLoad_Y[i] = -load[i][1]
            self.nodalLoad_Z[i] = -load[i][2]

# %% Other Functions

    def getNodalIndex(self,index):
        return self.FSI[index]

    def getNodalInitialPositions(self):
        return np.transpose(self.getPosition())

# %% Reads From the Lua File

    def read(self,path):

        file = open(path,'r')
        text = file.read().splitlines()
        file.close()

        for line in text:

            line = line.replace(' ','').replace('"','').replace("'",'')
            try: value = line.split('=')[1]
            except: continue

            if "Problem.id" in line: self.ID = value
            if "Problem.interface" in line: self.group = value
            if "Problem.maxFactor" in line: self.maxFactor = float(value)

# %% Other Functions

    def update(self,dt):

        self.mesh.remesh(False)
        if (self.ID == 'IncompNewtonNoT'): self.solver.precomputeMatrix()
        self.problem.copySolution(self.prevSolution)
        self.prevDisp = self.disp.copy()
        FluidSolver.update(self,dt)
        self.reload = False

    def save(self,_):
        self.problem.dump()

    def exit(self):

        print('======================================')
        self.problem.displayTimeStats()
        print('======================================')
        print('\nExit PFEM3D')

    def timeStats(self,time,dt):

        start = self.problem.getCurrentSimTime()
        print('[PFEM-1] : t1 = {:.5e} - dt = {:.3e}'.format(start,dt))
        print('[PFEM-2] : t2 = {:.5e} - factor = {:.0f}'.format(time,self.factor))
        print('----------------------------')