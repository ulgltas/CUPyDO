import cupydo.interfaces.Cupydo as cupy
from cupydo.testing import CTest,CTests
import numpy as np
import os

def test(meanFSIIt):

    name = [file for file in os.listdir() if('fluid' in file)]
    time = [float(file[8:-4]) for file in name]
    lastFile = name[np.argmax(time)]

    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber('General.Terminal',0)
    gmsh.open(lastFile)
    coord = gmsh.model.mesh.getNode(2)[0]

    tests = CTests()
    tests.add(CTest('Solid tip coordinate X',coord[0],0.320790,0.01,False))
    tests.add(CTest('Solid tip coordinate Y',coord[1],0.069850,0.01,False))
    tests.add(CTest('Mean number of ISI iterations',meanFSIIt,2.774387,0.01,False))
    tests.run()

# %% Input Parameters

def getFsiP():

    p = dict()
    path = os.path.abspath(os.path.dirname(__file__))
    p['cfdFile'] = path+'\input_pfem.lua'

    # Metafor and PFEM solvers
    
    p['fluidSolver'] = 'Pfem3D'
    p['solidSolver'] = 'Metafor'
    p['csdFile'] = 'input_meta'
    
    # FSI objects

    p['criterion'] = 'Displacements'
    p['interpolator'] = 'Matching'
    p['algorithm'] = 'IQN_ILS'
    
    # FSI parameters

    p['firstItTgtMat'] = False
    p['computation'] = 'direct'
    p['compType'] = 'unsteady'
    p['timeItTresh'] = 0
    p['dtSave'] = 0.01
    p['omega'] = 0.5
    p['nSteps'] = 0
    p['maxIt'] = 25
    p['tTot'] = 0.2
    p['tol'] = 1e-8
    p['dt'] = 1e-4
    p['nDim'] = 2
    
    return p

# %% Main Function

def main():

    param = getFsiP()
    cupydo = cupy.CUPyDO(param)
    cupydo.run()

    cupydo.algorithm.FluidSolver.save(cupydo.algorithm.timeIter)
    test(cupydo.algorithm.getMeanNbOfFSIIt())

if __name__=='__main__':
    main()