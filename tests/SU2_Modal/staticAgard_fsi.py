#!/usr/bin/env python3
# -*- coding: utf8 -*-

# CUPyDO configuration file
# Agard445 wing
# Adrien Crovato

def test(res, tol):
    import numpy as np
    from cupydo.testing import CTest, CTests
    # Read results from file
    with open("AerodynamicCoeff.ascii", 'rb') as f:
        lines = f.readlines()
    resultA = np.genfromtxt(lines[-1:], delimiter=None)
    with open("NodalDisplacement.dat", 'rb') as f:
        lines = f.readlines()
    resultS = np.genfromtxt(lines[-1:], delimiter=None)

    # Check convergence and results
    if (res > tol):
        print("\n\n" + "FSI residual = " + str(res) + ", FSI tolerance = " + str(tol))
        raise Exception("FSI algo failed to converge!")
    tests = CTests()
    tests.add(CTest('Lift coefficient', resultA[2], 0.053879, 1e-1, False))
    tests.add(CTest('Drag coefficient', resultA[3], 0.000032, 1e-1, False))
    tests.add(CTest('LE Displacement (16, z)', resultS[4], 0.0116, 1e-1, False))
    tests.add(CTest('TE Displacement (13808, z)', resultS[7], 0.0132, 1e-1, False))
    tests.run()

def getFsiP():
    """Fsi parameters"""
    import os
    filePath = os.path.abspath(os.path.dirname(__file__))
    p = {}
    # Solvers and config files
    p['fluidSolver'] = 'SU2'
    p['solidSolver'] = 'Modal'
    p['cfdFile'] = os.path.join(filePath,'staticAgard_fluid.cfg')
    p['csdFile'] = 'agard_solid'
    # FSI objects
    p['interpolator'] = 'RBF'
    p['interpType'] = 'conservative'
    p['algorithm'] = 'staticBGS'
    # FSI parameters
    p['regime'] = 'steady'
    p['computation'] = 'direct'
    p['nDim'] = 3
    p['dt'] = 0.
    p['tTot'] = 0.
    
    p['dtSave'] = 0
    p['tol'] = 5e-3
    p['maxIt'] = 6
    p['omega'] = 1.0
    p['rbfRadius'] = 1.
    p['nodalLoadsType'] = 'force'
    # Coupling Type

    p['mechanical'] = True
    p['thermal'] = False
    return p

def main():
    import cupydo.interfaces.Cupydo as cupy
    p = getFsiP() # get parameters
    cupydo = cupy.CUPyDO(p) # create fsi driver
    cupydo.run() # run fsi process
    test(cupydo.algorithm.criterion.epsilon, p['tol']) # check the results
    
    # eof
    print('')

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':
    main()
