#!/usr/bin/env python
# -*- coding: utf-8 -*-

# CUPyDO configuration file
# Agard445 wing
# Adrien Crovato & Mariano Sánchez Martínez

def test(cupydo, tol):
    res = cupydo.algorithm.errValue
    import numpy as np
    from cupydo.testing import *
    # Read results from data
    cl = cupydo.algorithm.FluidSolver.coreSolver.getCl()
    cd = cupydo.algorithm.FluidSolver.coreSolver.getCd()
    with open("db_Field(TZ,RE)_GROUP_ID_121.ascii", 'rb') as f:
        lines = f.readlines()
    resultS1 = np.genfromtxt(lines[-1:], delimiter=None)
    with open("db_Field(TZ,RE)_GROUP_ID_122.ascii", 'rb') as f:
        lines = f.readlines()
    resultS2 = np.genfromtxt(lines[-1:], delimiter=None)

    # Check convergence and results
    if (res > tol):
        print "\n\n" + "FSI residual = " + str(res) + ", FSI tolerance = " + str(tol)
        raise Exception(ccolors.ANSI_RED + "FSI algo failed to converge!" + ccolors.ANSI_RESET)
    tests = CTests()
    tests.add(CTest('Lift coefficient', cl, 0.0460, 1e-2, True)) # abs. tol
    tests.add(CTest('Drag coefficient', cd, 0.00080, 1e-4, True))
    tests.add(CTest('LE vertical displacement', resultS2[2], 0.009, 5e-3, True))
    tests.add(CTest('TE vertical displacement', resultS1[2], 0.010, 5e-3, True))
    tests.run()

def getFsiP():
    """Fsi parameters"""
    import os
    fileName = os.path.splitext(os.path.basename(__file__))[0]
    p = {}
    # Solvers and config files
    p['fluidSolver'] = 'VLM'
    p['solidSolver'] = 'Metafor'
    p['cfdFile'] = fileName[:-3] + 'fluid'
    p['csdFile'] = fileName[:-3] + 'solid'
    # FSI objects
    p['interpolator'] = 'RBF'
    p['criterion'] = 'Displacements'
    p['algorithm'] = 'StaticBGS'
    # FSI parameters
    p['compType'] = 'steady'
    p['nDim'] = 3
    p['dt'] = 0.1
    p['tTot'] = 0.1
    p['timeItTresh'] = -1
    p['tol'] = 1e-4
    p['maxIt'] = 50
    p['omega'] = 1.0
    p['rbfRadius'] = .5
    return p

def main():
    import cupydo.interfaces.Cupydo as cupy
    import staticAgard_fluid
    p = getFsiP() # get parameters
    cupydo = cupy.CUPyDO(p) # create fsi driver
    cupydo.run() # run fsi process
    cupydo.algorithm.FluidSolver.coreSolver.save()
    test(cupydo, p['tol']) # check the results
    
    # eof
    print ''

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':
    main()
