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
        path = param['cfdFile']+'.lua'
        input = read(path)

        # Problem type availables

        IC = ['IncompNewtonNoT','Bingham','Boussinesq','BoussinesqPC','Conduction','ConductionPC']
        WC = ['WCompNewtonNoT','WCBoussinesq','WCBoussinesqPC']

        # Problem initialization

        if not input['useCupydo']: raise Exception('UseCupydo must be True')
        elif input['id'] in IC: self.problem = w.ProbIncompNewton(path)
        elif input['id'] in WC: self.problem = w.ProbWCompNewton(path)
        else: raise Exception('Unknown problem type '+input['id'])

        # Gets some important objects and variables

        self.solver = self.problem.getSolver()
        self.mesh = self.problem.getMesh()

        # Initializes the load and mesh copy

        self.meshSave = w.Mesh()
        self.FSID = w.VectorInt()
        self.load = w.VectorArrayDouble3()

        # Gets Physical nodes at the FS interface

        self.mesh.getNodesIndexTag('FSInterface',self.FSID)
        self.nPhysicalNodes = self.FSID.size()
        self.nNodes = self.FSID.size()

        # Initializes the fluid solver

        self.problem.updateTime(-param['dt'])
        FluidSolver.__init__(self)
        self.initRealTimeData()
        
# %% Calculates One Increment From t1 to t2

    def run(self,t1,t2):
        
        if self.reload:

            self.problem.loadMesh(self.meshSave,t1)
            self.applyNodalDisplacements(*self.disp,0,0,0,0,0)

        # Initialization of the simulation

        self.mesh.remesh(self.problem.isOutputVerbose())
        time = self.problem.getCurrentSimTime()
        dt = self.solver.getTimeStep()
        progress = True
        print()

        # Solves without remeshing until t2

        while progress:
            if (dt+time-t2)/dt > -0.2:
                    
                self.solver.setTimeStep(t2-time)
                ok = self.solver.solveOneTimeStep()
                if(ok): progress = False

            else: ok = self.solver.solveOneTimeStep()
            time = self.problem.getCurrentSimTime()
            dt = self.solver.getTimeStep()
            self.solver.computeNextDT()

            # Prints the current state

            print(ok,': t = {:.6f} - dt = {:.2e}'.format(time,dt))
            if(ok): self.mesh.remesh(self.problem.isOutputVerbose())
            dt = self.solver.getTimeStep()

        # Computes nodal fluid loads

        self.setCurrentState()
        self.reload = True

# %% Get and Set Nodal Values

    def getNodalInitialPositions(self):

        coord = np.zeros((self.nPhysicalNodes,3))
        for i in range(self.nPhysicalNodes):

            idx = self.FSID[i]
            node = self.mesh.getNode(idx)
            coord[i] = [node.getCoordinate(j) for j in range(3)]

        return coord.T

    # Sets current node states

    def setCurrentState(self):
        
        self.mesh.getNodesIndexTag('FSInterface',self.FSID)
        self.solver.computeLoads(self.load)
        for i in range(self.nPhysicalNodes):

            idx = self.FSID[i]
            self.nodalLoad_X[i] = self.load[idx][0]
            self.nodalLoad_Y[i] = self.load[idx][1]
            self.nodalLoad_Z[i] = self.load[idx][2]

    # Returns the index of the index-th interface node

    def getNodalIndex(self,index):
        return self.FSID[index]

    # Prescribes Nodal Positions from Solid Solver

    def applyNodalDisplacements(self,dx,dy,dz,dx_nM1,dy_nM1,dz_nM1,haloNodesDisplacements,time):

        self.disp = [dx,dy,dz]
        disp = np.add(self.disp,self.disp0).T
        self.mesh.getNodesIndexTag('FSInterface',self.FSID)
        position = w.VectorArrayDouble3(self.nPhysicalNodes)
        for i in range(self.nPhysicalNodes): position[i] = disp[i]
        self.mesh.setNodesPosition(position,self.FSID)

# %% Update and Save Results

    def update(self,dt):

        self.reload = False
        self.problem.copyMesh(self.meshSave)
        FluidSolver.update(self,dt)

    def save(self,nt):
        self.problem.writeExtractors()

    def initRealTimeData(self):

        self.problem.dump()
        self.problem.displayParams()
        self.disp0 = self.getNodalInitialPositions()
        self.reload = False
        self.time = -np.inf

# %% Exits the Solver

    def exit(self):

        print('======================================')
        self.problem.displayTimeStats()
        print('======================================')
        print('\n***************************** Exit PFEM3D *****************************')
        return