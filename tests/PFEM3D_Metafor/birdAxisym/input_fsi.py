import cupydo.interfaces.Cupydo as cupy
from cupydo.testing import CTest,CTests
import numpy as np
import gmsh
import os

def test(meanFSIIt):

    name = [file for file in os.listdir() if('solid' in file)]
    time = [float(file[8:-4]) for file in name]
    lastFile = name[np.argmax(time)]

    if not gmsh.isInitialized(): gmsh.initialize()
    gmsh.option.setNumber('General.Terminal',0)
    gmsh.open(lastFile)
    coord = gmsh.model.mesh.getNode(4)[0]
    if gmsh.isInitialized(): gmsh.finalize()

    tests = CTests()
    tests.add(CTest('Middle bar coordinate X', coord[0], 0, 1e-9, True))
    tests.add(CTest('Middle bar coordinate Y', coord[1], -0.000894, 0.05, False))
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
    p['interpType'] = 'conservative'
    p['algorithm'] = 'aitkenBGS'
    
    # FSI parameters

    p['firstItTgtMat'] = False
    p['computation'] = 'direct'
    p['regime'] = 'unsteady'
    p['dtSave'] = 0
    p['omega'] = 0.5
    p['maxIt'] = 25
    p['dt'] = 1.5e-6
    p['tTot'] = 1e-4
    p['criterion'] = 'relative'
    p['nDim'] = 2

    # Coupling Type

    p['mechanical'] = True
    p['mechanicalTol'] = 1e-6
    p['thermal'] = False
    return p

# Main Function

def main():

    param = getFsiP()
    cupydo = cupy.CUPyDO(param)
    cupydo.run()

    test(cupydo.algorithm.getMeanNbOfFSIIt())

if __name__=='__main__':
    main()
