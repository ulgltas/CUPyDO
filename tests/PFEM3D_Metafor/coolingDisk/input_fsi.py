import cupydo.interfaces.Cupydo as cupy
from cupydo.testing import CTest,CTests
import numpy as np
import gmsh
import os

def test(meanFSIIt):

    output = np.loadtxt('output.txt',delimiter=',')
    temperature = output[-1][-1]

    name = [file for file in os.listdir() if('solid' in file)]
    time = [float(file[8:-4]) for file in name]
    lastFile = name[np.argmax(time)]
    tag = 73

    if not gmsh.isInitialized(): gmsh.initialize()
    gmsh.option.setNumber('General.Terminal',0)
    gmsh.open(lastFile)
    coord = gmsh.model.mesh.getNode(tag)[0]
    if gmsh.isInitialized(): gmsh.finalize()

    tests = CTests()
    tests.add(CTest('Center ball coordinate X', coord[0], 0.138, 1e-3, False))
    tests.add(CTest('Center ball coordinate Y', coord[1], 0.172, 0.01, False))
    tests.add(CTest('Center ball temperature', temperature, 168.2, 0.005, False))
    tests.add(CTest('Mean number of ISI iterations', meanFSIIt, 9, 1, True))
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

    p['interpolator'] = 'RBF'
    p['interpType'] = 'consistent'
    p['chtTransferMethod'] = 'FFTB'
    p['algorithm'] = 'IQN_MVJ'
    p['rbfRadius'] = 0.1
    
    # FSI parameters

    p['firstItTgtMat'] = False
    p['computation'] = 'direct'
    p['regime'] = 'unsteady'
    p['dtSave'] = 1e-1
    p['omega'] = 0.5
    p['maxIt'] = 25
    p['tol'] = 1e-4
    p['dt'] = 1e-2
    p['tTot'] = 1
    p['nDim'] = 2

    # Coupling Type

    p['mechanical'] = True
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
