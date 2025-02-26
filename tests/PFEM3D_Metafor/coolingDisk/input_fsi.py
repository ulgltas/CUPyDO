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
    tests.add(CTest('Center ball coordinate X', coord[0], 0.15, 0.1, False))
    tests.add(CTest('Center ball coordinate Y', coord[1], 0.16, 0.1, False))
    tests.add(CTest('Center ball temperature', temperature, 170, 0.01, False))
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
    p['dt'] = 5e-3
    p['tTot'] = 1
    p['criterion'] = 'relative'
    p['nDim'] = 2
    p['qrFilter'] = None

    # PETSC parameters

    p['interpMaxIt'] = 100
    p['interpTol'] = 1e-8

    # Coupling Type

    p['mechanical'] = True
    p['mechanicalTol'] = 1e-6
    p['thermal'] = True
    p['thermalTol'] = 1e-6
    return p

# Main Function

def main():

    param = getFsiP()
    cupydo = cupy.CUPyDO(param)
    cupydo.run()

    test(cupydo.algorithm.getMeanNbOfFSIIt())

if __name__=='__main__':
    main()
