#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# CUPyDO configuration file
# Agard445 wing
# Adrien Crovato

def test(res, tol):
    import numpy as np
    from cupydo.testing import CTest, CTests
    # Read results from file
    with open("DartHistory.dat", 'rb') as f:
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
    tests.add(CTest('Lift coefficient', resultA[2], 0.0537, 1e-2, True)) # abs. tol
    tests.add(CTest('Drag coefficient', resultA[3], 0.00045, 1e-4, True))
    tests.add(CTest('LE Displacement (16, z)', resultS[4], 0.0116, 1e-1, False)) # rel. tol. of 10%
    tests.add(CTest('TE Displacement (13808, z)', resultS[7], 0.0132, 1e-1, False)) # rel. tol. of 10%
    tests.run()

def getFsiP():
    """Fsi parameters"""
    import os
    fileName = os.path.splitext(os.path.basename(__file__))[0]
    p = {}
    
    # Solvers and config files

    p['fluidSolver'] = 'DART'
    p['solidSolver'] = 'Modal'
    p['cfdFile'] = fileName[:-3] + 'fluid'
    p['csdFile'] = fileName[:-3] + 'solid'

    # FSI objects

    p['interpolator'] = 'RBF'
    p['interpType'] = 'conservative'
    p['algorithm'] = 'IQN_ILS'

    # FSI parameters

    p['regime'] = 'steady'
    p['computation'] = 'direct'
    p['criterion'] = 'relative'
    p['nDim'] = 3
    p['dt'] = 0.0
    p['tTot'] = 0.0
    p['dtSave'] = 0
    p['maxIt'] = 50
    p['omega'] = 1.0
    p['nSteps'] = 0
    p['firstItTgtMat'] = False
    p['rbfRadius'] = 1.0
    p['qrFilter'] = 'Haelterman'

    # Coupling Type

    p['mechanical'] = True
    p['mechanicalTol'] = 1e-5
    p['thermal'] = False
    return p

def main():
    import cupydo.interfaces.Cupydo as cupy
    p = getFsiP() # get parameters
    cupydo = cupy.CUPyDO(p) # create fsi driver
    cupydo.run() # run fsi process
    test(cupydo.algorithm.criterion.epsilon, p['mechanicalTol']) # check the results
    
    # eof
    print('')

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':
    main()