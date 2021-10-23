import cupydo.interfaces.Cupydo as cupy
from cupydo.testing import CTest,CTests
import numpy as np
import os

def test(mean):

    name = [file for file in os.listdir() if('fluid' in file)]
    time = [float(file[8:-4]) for file in name]
    lastFile = name[np.argmax(time)]

    import gmsh
    gmsh.initialize()
    gmsh.open(lastFile)
    coord,_ = gmsh.model.mesh.getNode(2)

    tests = CTests()
    tests.add(CTest('Mean nb of FSI iterations',mean,1.586139,0.1,True))
    tests.add(CTest('Solid tip coordinate X',coord[0],0.311551,1e-3,True))
    tests.add(CTest('Solid tip coordinate Y',coord[1],0.070568,1e-4,True))
    tests.run()

# %% Input Parameters

def getFsiP():

    param = {}
    path = os.path.abspath(os.path.dirname(__file__))
    param['cfdFile'] = path+'\\'

    # Metafor and PFEM solvers
    
    param['fluidSolver'] = 'Pfem3D'
    param['solidSolver'] = 'Metafor'
    param['cfdFile'] += 'input_pfem'
    param['csdFile'] = 'input_meta'
    
    # FSI objects

    param['criterion'] = 'Displacements'
    param['interpolator'] = 'Matching'
    param['algorithm'] = 'AitkenBGS'
    
    # FSI parameters

    param['computation'] = 'direct'
    param['compType'] = 'unsteady'
    param['timeItTresh'] = 0
    param['omega'] = 0.5
    param['tTot'] = 0.3
    param['tol'] = 1e-8
    param['dt'] = 1e-5
    param['maxIt'] = 5
    param['nDim'] = 2
    
    return param

# %% Main Function

def main():

    param = getFsiP()
    cupydo = cupy.CUPyDO(param)
    cupydo.run()

    mean = cupydo.algorithm.getMeanNbOfFSIIt()
    test(mean)

if __name__=='__main__':
    main()