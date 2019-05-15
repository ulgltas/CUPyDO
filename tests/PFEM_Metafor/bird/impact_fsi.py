
#! /usr/bin/env python
# -*- coding: utf-8 -*-
# original name: birdImpact_deformable_panel_alu_Mtf_Pfem_fsi

def test(res, tol, it):
    import numpy as np
    from cupydo.testing import *
    # Check convergence and results
    if (algorithm.errValue > p['tollFSI']):
        print "\n\n" + "FSI residual = " + \
            str(algorithm.errValue) + ", FSI tolerance = " + str(p['tollFSI'])
        raise Exception(ccolors.ANSI_RED +
                        "FSI algo failed to converge!" + ccolors.ANSI_RESET)

    # Read results from file
    with open("db_Field(TY,RE)_GROUP_ID_17.ascii", 'rb') as f:
        lines = f.readlines()
    result_1 = np.genfromtxt(lines[-1:], delimiter=None)

    tests = CTests()
    tests.add(CTest('Mean nb of FSI iterations', it, 2, 1, True))
    tests.add(CTest('Y-displacement panel center', result_1[2], -0.001462, 1e-2, False))
    tests.run()

def getFsiP():
    """Fsi parameters"""
    import os
    fileName = os.path.splitext(os.path.basename(__file__))[0]
    p = {}
    # For this case
    U0 = 100
    N = 10
    R = 0.01
    d = 2.5*R/N
    # Solvers and config files
    p['fluidSolver'] = 'Pfem'
    p['solidSolver'] = 'Metafor'
    p['cfdFile'] = fileName[:-3] + 'fluid'
    p['csdFile'] = fileName[:-3] + 'solid'
    # FSI objects
    p['interpolator'] = 'Matching'
    p['criterion'] = 'Displacements'
    p['algorithm'] = 'AitkenBGS'
    # FSI parameters
    p['compType'] = 'unsteady'
    p['nDim'] = 2
    p['dt'] = 1.5e-6
    p['tTot'] = 1e-4 # 40*((4*R)/U0 + d/U0)
    p['timeItTresh'] = 0
    p['tol'] = 1e-6
    p['maxIt'] = 20
    p['omega'] = 0.5
    p['nBnd'] = 13
    p['saveFreqPFEM'] = 10
    p['mtfSaveAllFacs'] = False
    return p

def main():
    import cupydo.interfaces.CUPYDO as cupy
    p = getFsiP() # get parameters
    cupydo = cupy.Cupydo(p) # create fsi driver
    cupydo.run() # run fsi process
    test(cupydo.algorithm.errValue, p['tol'], cupydo.algorithm.getMeanNbOfFSIIt()) # check the results
    
    # eof
    print ''

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':
    main()
