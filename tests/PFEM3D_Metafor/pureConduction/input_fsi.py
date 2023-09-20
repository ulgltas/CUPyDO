import cupydo.interfaces.Cupydo as cupy
from cupydo.testing import CTest,CTests
import numpy as np
import os

def test(meanFSIIt):

    output = np.loadtxt('output.txt',delimiter=',')
    temperature = output[-1][-1]

    tests = CTests()
    tests.add(CTest('Middle bar temperature', temperature, 311.75, 0.005, False))
    tests.add(CTest('Mean number of ISI iterations', meanFSIIt, 5, 1, True))
    tests.run()

# Input Parameters

def getFsiP():

    p = dict()
    path = os.path.abspath(os.path.dirname(__file__))
    p['cfdFile'] = path+'/input_pfem.lua'

    # Metafor and PFEM solvers
    
    p['fluidSolver'] = 'Pfem3D'
    p['solidSolver'] = 'Metafor'
    p['csdFile'] = 'input_meta'
    
    # FSI objects

    p['interpolator'] = 'matching'
    p['interpType'] = 'consistent'
    p['chtTransferMethod'] = 'FFTB'
    p['algorithm'] = 'IQN_ILS'
    p['nSteps'] = 5
    
    # FSI parameters

    p['firstItTgtMat'] = False
    p['computation'] = 'direct'
    p['compType'] = 'unsteady'
    p['timeItTresh'] = 0
    p['dtSave'] = 0.1
    p['omega'] = 0.5
    p['maxIt'] = 25
    p['tol'] = 1e-6
    p['dt'] = 0.1
    p['tTot'] = 20
    p['nDim'] = 2

    # Coupling Type

    p['mechanical'] = False
    p['thermal'] = True
    return p

# Main Function

def main():

    param = getFsiP()
    cupydo = cupy.CUPyDO(param)
    cupydo.run()

    test(cupydo.algorithm.getMeanNbOfFSIIt())

if __name__=='__main__':
    main()
