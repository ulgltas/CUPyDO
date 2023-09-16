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
    with open("db_Field(TZ,RE)_GROUP_ID_121.ascii", 'rb') as f:
        lines = f.readlines()
    resultS1 = np.genfromtxt(lines[-1:], delimiter=None)
    with open("db_Field(TZ,RE)_GROUP_ID_122.ascii", 'rb') as f:
        lines = f.readlines()
    resultS2 = np.genfromtxt(lines[-1:], delimiter=None)

    # Check convergence and results
    if (res > tol):
        print("\n\n" + "FSI residual = " + str(res) + ", FSI tolerance = " + str(tol))
        raise Exception("FSI algo failed to converge!")
    tests = CTests()
    tests.add(CTest('Lift coefficient', resultA[2], 0.0537, 1e-2, True)) # abs. tol
    tests.add(CTest('Drag coefficient', resultA[3], 0.00045, 1e-4, True))
    tests.add(CTest('LE vertical displacement', resultS2[-1], 0.0116, 5e-3, True))
    tests.add(CTest('TE vertical displacement', resultS1[-1], 0.0132, 5e-3, True))
    tests.run()

def getFsiP():
    """Fsi parameters"""
    import os
    fileName = os.path.splitext(os.path.basename(__file__))[0]
    p = {}
    # Solvers and config files
    p['fluidSolver'] = 'DART'
    p['solidSolver'] = 'Metafor'
    p['cfdFile'] = fileName[:-3] + 'fluid'
    p['csdFile'] = fileName[:-3] + 'solid'
    # FSI objects
    p['interpolator'] = 'RBF'
    p['interpType'] = 'conservative'
    p['algorithm'] = 'staticBGS'
    # FSI parameters
    p['compType'] = 'steady'
    p['computation'] = 'direct'
    p['nDim'] = 3
    p['dt'] = 0.1
    p['tTot'] = 0.1
    
    p['dtSave'] = 0
    p['tol'] = 1e-5
    p['maxIt'] = 50
    p['omega'] = 1.0
    p['rbfRadius'] = .3
    # Coupling Type

    p['mechanical'] = True
    p['thermal'] = False
    return p

def main():
    import cupydo.interfaces.Cupydo as cupy
    p = getFsiP() # get parameters
    cupydo = cupy.CUPyDO(p) # create fsi driver
    cupydo.run() # run fsi process
    test(cupydo.algorithm.errValue, p['tol']) # check the results
    
    # eof
    print('')

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':
    main()
