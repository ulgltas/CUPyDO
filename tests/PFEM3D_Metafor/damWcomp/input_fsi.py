import cupydo.interfaces.Cupydo as cupy
from cupydo.testing import CTest,CTests
import numpy as np
import os

def test():

    name = [file for file in os.listdir() if('fluid' in file)]
    time = [float(file[8:-4]) for file in name]
    lastFile = name[np.argmax(time)]

    import gmsh
    gmsh.initialize()
    gmsh.open(lastFile)
    coord = gmsh.model.mesh.getNode(2)[0]

    tests = CTests()
    tests.add(CTest('Solid tip coordinate X',coord[0],0.305965,0.1,False))
    tests.add(CTest('Solid tip coordinate Y',coord[1],0.070794,0.1,False))
    tests.run()

# %% Input Parameters

def getFsiP():

    param = {}
    path = os.path.abspath(os.path.dirname(__file__))
    param['cfdFile'] = path+'\input_pfem.lua'

    # Metafor and PFEM solvers
    
    param['fluidSolver'] = 'Pfem3D'
    param['solidSolver'] = 'Metafor'
    param['csdFile'] = 'input_meta'
    
    # FSI objects

    param['criterion'] = 'Displacements'
    param['interpolator'] = 'Matching'
    param['algorithm'] = 'IQN_ILS'
    
    # FSI parameters

    param['firstItTgtMat'] = False
    param['computation'] = 'direct'
    param['compType'] = 'unsteady'
    param['timeItTresh'] = 0
    param['dtWrite'] = 0.01
    param['omega'] = 0.5
    param['nSteps'] = 0
    param['maxIt'] = 25
    param['tTot'] = 0.3
    param['tol'] = 1e-8
    param['dt'] = 1e-4
    param['nDim'] = 2
    
    return param

# %% Main Function

def main():

    param = getFsiP()
    cupydo = cupy.CUPyDO(param)
    cupydo.run()
    test()

if __name__=='__main__':
    main()