import cupydo.interfaces.Cupydo as cupy
from cupydo.testing import CTest,CTests
import os

def test(mean):
    
    tests = CTests()
    tests.add(CTest('Mean nb of FSI iterations',mean,1.6004,1e-5,True))
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

    param['interpolator'] = 'Matching'
    param['criterion'] = 'Displacements'
    param['algorithm'] = 'AitkenBGS'
    
    # FSI parameters

    param['compType'] = 'unsteady'
    param['computation'] = 'direct'
    param['timeItTresh'] = 0
    param['nDim'] = 2
    param['dt'] = 1e-4
    param['tTot'] = 0.3
    param['tol'] = 1e-6
    param['maxIt'] = 25
    param['omega'] = 0.5
    
    return param

def main():

    param = getFsiP()
    cupydo = cupy.CUPyDO(param)
    cupydo.run()

    mean = cupydo.algorithm.getMeanNbOfFSIIt()
    test(mean)

if __name__=='__main__':
    main()