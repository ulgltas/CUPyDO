#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# CUPyDO configuration file
# Diamond airfoil
# Adrien Crovato

def test(res, tol):
    import numpy as np
    from cupydo.testing import CTest, CTests, ccolors
    # Read results from file
    with open("DartHistory.dat", 'rb') as f:
        lines = f.readlines()
    resultA = np.genfromtxt(lines[-1:], delimiter=None)
    with open("db_Field(TY,RE)_GROUP_ID_101.ascii", 'rb') as f:
        lines = f.readlines()
    resultS = np.genfromtxt(lines[-1:], delimiter=None)

    # Check convergence and results
    if (res > tol):
        print("\n\n" + "FSI residual = " + str(res) + ", FSI tolerance = " + str(tol))
        raise Exception(ccolors.ANSI_RED + "FSI algo failed to converge!" + ccolors.ANSI_RESET)
    tests = CTests()
    tests.add(CTest('Lift coefficient', resultA[2], 0.22, 1e-1, False)) # rel. tol. of 10%, dummy value
    tests.add(CTest('TE. vertical displacement', resultS[2], 0.034, 1e-1, False)) # rel. tol. of 10%, dummy value
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
    p['criterion'] = 'Displacements'
    p['algorithm'] = 'IQN_ILS'
    # FSI parameters
    p['compType'] = 'steady'
    p['computation'] = 'direct'
    p['nDim'] = 2
    p['dt'] = 0.1
    p['tTot'] = 0.1
    p['timeItTresh'] = -1
    p['dtSave'] = 0
    p['tol'] = 1e-5
    p['maxIt'] = 50
    p['omega'] = 1.0
    p['nSteps'] = 0
    p['firstItTgtMat'] = False
    p['rbfRadius'] = 1.
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
