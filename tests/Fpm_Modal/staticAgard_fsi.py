#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# CUPyDO configuration file
# Agard445 wing
# Adrien Crovato & Axel Dechamps

def test(res, tol):
    import numpy as np
    from cupydo.testing import CTest, CTests, ccolors
    # Read results from file
    with open("FlowHistory.dat", 'rb') as f:
        lines = f.readlines()
    resultA = np.genfromtxt(lines[-1:], delimiter=None)
    with open("NodalDisplacement.dat", 'rb') as f:
        lines = f.readlines()
    resultS = np.genfromtxt(lines[-1:], delimiter=None)

    # Check convergence and results
    if (res > tol):
        print ("\n\n" + "FSI residual = " + str(res) + ", FSI tolerance = " + str(tol))
        raise Exception(ccolors.ANSI_RED + "FSI algo failed to converge!" + ccolors.ANSI_RESET)
    tests = CTests()
    tests.add(CTest('Lift coefficient', resultA[2], 0.0459, 1e-2, True)) # Flow M=0: 0.0447
    tests.add(CTest('Drag coefficient', resultA[3], 0.00025, 1e-3, True)) # Flow M=0: 0.0004
    tests.add(CTest('LE Displacement (16, z)', resultS[4], 0.0094, 1e-1, False)) # rel. tol. of 10%
    tests.add(CTest('TE Displacement (13808, z)', resultS[7], 0.0098, 1e-1, False)) # rel. tol. of 10%
    tests.run()

def getFsiP():
    """Fsi parameters"""
    import os
    fileName = os.path.splitext(os.path.basename(__file__))[0]
    p = {}
    # Solvers and config files
    p['computation'] = 'Direct'             # type of computation available: Adjoint or Direct
    p['fluidSolver'] = 'Fpm'                # fluid solvers available: SU2, Pfem, Flow, Fpm 
    p['solidSolver'] = 'Modal'              # solid solvers available: Metafor, RBMI Modal, GetDP 
    p['cfdFile'] = fileName[:-3] + 'fluid'  # path to fluid cfg file 
    p['csdFile'] = fileName[:-3] + 'solid'  # path to solid cfg file 
    # FSI objects
    p['interpolator'] = 'RBF'               # interpolator type available: Matching, RBF, TPS 
    p['criterion'] = 'Displacements'        # convergence criterion available: Displacements 
    p['algorithm'] = 'IQN_ILS'              # FSI algorithms available: Explicit, StaticBGS, AitkenBGS, IQN_ILS 
    # FSI parameters
    p['compType'] = 'steady'                # steady or unsteady 
    p['nDim'] = 3                           # dimension: 2 or 3 
    p['dt'] = 0.0                           # time steps (if unsteady) 
    p['tTot'] = 0.0                         # total time
    p['timeItTresh'] = -1                   # ?
    p['tol'] = 1e-5                         # tolerance on displacements (required by BGS)
    p['maxIt'] = 50                         # maximum number of iterations (required by BGS)
    p['omega'] = 1.0                        # relaxation parameter (required by BGS)
    p['nSteps'] = 0                         # number of time steps to keep (required by IQN_ILS)
    p['firstItTgtMat'] = False              # compute the tangent matrix based on first iteration: true or false (required by IQN_ILS)
    p['rbfRadius'] = 1.                     # radius of interpolation (required by RBF)
    return p

def main():
    import cupydo.interfaces.Cupydo as cupy
    p = getFsiP() # get parameters
    cupydo = cupy.CUPyDO(p) # create fsi driver
    cupydo.run() # run fsi process
    test(cupydo.algorithm.errValue, p['tol']) # check the results
    
    # eof
    print ('')

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':
    main()