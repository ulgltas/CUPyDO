import cupydo.interfaces.Cupydo as cupy
from cupydo.testing import CTest,CTests
import numpy as np
import gmsh
import os

def test(meanFSIIt):

    name = [file for file in os.listdir() if('solid' in file)]
    time = [float(file[8:-4]) for file in name]
    lastFile = name[np.argmax(time)]
    tag = 2

    if not gmsh.isInitialized(): gmsh.initialize()
    gmsh.option.setNumber('General.Terminal',0)
    gmsh.open(lastFile)
    coord = gmsh.model.mesh.getNode(tag)[0]
    if gmsh.isInitialized(): gmsh.finalize()

    tests = CTests()
    tests.add(CTest('Solid tip coordinate X', coord[0], 0.333640, 0.05, False))
    tests.add(CTest('Solid tip coordinate Y', coord[1], 0.071319, 0.05, False))
    tests.add(CTest('Mean number of ISI iterations', meanFSIIt, 4, 1, True))
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
    p['algorithm'] = 'aitkenBGS'
    
    # FSI parameters

    p['firstItTgtMat'] = False
    p['computation'] = 'direct'
    p['regime'] = 'unsteady'
    
    p['omega'] = 0.5
    p['dtSave'] = 0.01
    p['maxIt'] = 20
    p['tTot'] = 0.35
    p['tol'] = 1e-6
    p['dt'] = 0.001
    p['nDim'] = 2
    
    # Coupling Type

    p['mechanical'] = True
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