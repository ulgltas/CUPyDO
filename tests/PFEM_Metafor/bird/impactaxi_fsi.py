#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# original name: birdImpact_deformable_panel_alu_Mtf_Pfem_Axisym_fsi.py

def test(res, tol, it):
    import numpy as np
    from cupydo.testing import CTest, CTests
    # Check convergence and results
    if (res > tol):
        print("\n\n" + "FSI residual = " + str(res) + ", FSI tolerance = " + str(tol))
        raise Exception("FSI algo failed to converge!")
    
    # Read results from file
    with open("db_Field(TY,RE)_GROUP_ID_18.ascii", 'rb') as f:
        lines = f.readlines()
    result_1 = np.genfromtxt(lines[-1:], delimiter=None)
    
    tests = CTests()
    tests.add(CTest('Mean nb of FSI iterations', it, 5, 1, True))
    tests.add(CTest('Y-displacement panel center', result_1[-1], -0.000894, 0.05, False))
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

    p['interpolator'] = 'matching'
    p['interpType'] = 'conservative'
    p['algorithm'] = 'aitkenBGS'

    # FSI parameters

    p['regime'] = 'unsteady'
    p['computation'] = 'direct'
    p['criterion'] = 'relative'
    p['nDim'] = 2
    p['dt'] = 1.5e-6
    p['tTot'] = 1e-4 # 40*((4*R)/U0 + d/U0)
    
    p['dtSave'] = 0
    p['maxIt'] = 20
    p['omega'] = 0.5

    # Coupling Type

    p['mechanical'] = True
    p['mechanicalTol'] = 1e-6
    p['thermal'] = False
    return p

def main():
    import cupydo.interfaces.Cupydo as cupy
    p = getFsiP() # get parameters
    cupydo = cupy.CUPyDO(p) # create fsi driver
    cupydo.run() # run fsi process
    test(cupydo.algorithm.criterion.epsilon, p['mechanicalTol'], cupydo.algorithm.getMeanNbOfFSIIt()) # check the results
    
    # eof
    print('')

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':
    main()