#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# CUPyDO configuration file
# Naca0012 airfoil
# Adrien Crovato

def test(res, tol):
    import numpy as np
    from cupydo.testing import CTest, CTests
    # Flow constant (defined in case_fluid)
    dynP = 0.5*100 # dynamic pressure
    alpha0 = np.radians(3) # initial angle of attack
    cRef = 1

    # RBM constant (defined in case_solid)
    k = 250 # vertical stiffness
    kappa = 25 # rotational stiffness

    # Airfoil slopes (measured from DART)
    cl_alpha = 6.8
    cd_alpha = 0.085
    cm_alpha = 0.27 # must be measured from flexural axis posit. hard to calibrate but of crucial importance!

    # Pitch-plunge system of equations
    A = np.array([ [k, -dynP*cRef*(cl_alpha*np.cos(alpha0) - cd_alpha*np.sin(alpha0))],
                   [0, kappa - dynP*cRef*cRef*cm_alpha] ])
    b = np.array([0, kappa*alpha0])
    x = np.linalg.solve(A, b)

    # Display the solution
    print("Ref. lift coefficient: " + str(cl_alpha*x[1]))
    print("Ref. vertical displacement: " + str(x[0]))
    print("Ref. new angle of attack: " + str(np.degrees(x[1])))
    print("Ref. rotational displacement : " + str(np.degrees(x[1]-alpha0)))

    # Read results from file
    resultS = np.genfromtxt("NativeHistory.dat", delimiter=None, skip_header=1)
    with open("DartHistory.dat", 'rb') as f:
        lines = f.readlines()
    resultA = np.genfromtxt(lines[-1:], delimiter=None)

    # Check convergence and results
    if (res > tol):
        print("\n\n" + "FSI residual = " + str(res) + ", FSI tolerance = " + str(tol))
        raise Exception("FSI algo failed to converge!")
    tests = CTests()
    tests.add(CTest('Lift coefficient', resultA[2], cl_alpha*x[1], 2*1e-1, False)) # rel. tol. of 10%
    tests.add(CTest('Vertical displacement', -resultS[2]/cRef, x[0]/cRef, 2*1e-2, True)) # abs. tol. of 1% of chord
    tests.add(CTest('Rotational displacement', np.degrees(resultS[3]), np.degrees(x[1]-alpha0), 2*5e-1, True)) # abs. tol. of .5°
    tests.run()

def getFsiP():
    """Fsi parameters"""
    import os
    filePath = os.path.abspath(os.path.dirname(__file__))
    fileName = os.path.splitext(os.path.basename(__file__))[0]
    p = {}
    # Solvers and config files
    p['fluidSolver'] = 'DART'
    p['solidSolver'] = 'RBMI'
    p['cfdFile'] = fileName[:-3] + 'fluid'
    p['csdFile'] = os.path.join(filePath, fileName[:-3] + 'solid.cfg')
    # FSI objects
    p['interpolator'] = 'matching'
    p['criterion'] = 'displacement'
    p['algorithm'] = 'IQN_ILS'
    # FSI parameters
    p['compType'] = 'Steady'
    p['computation'] = 'direct'
    p['nDim'] = 2
    p['dt'] = 0.0
    p['tTot'] = 0.0
    
    p['dtSave'] = 0
    p['tol'] = 1e-4
    p['maxIt'] = 50
    p['omega'] = 1.0
    p['nSteps'] = 0
    p['firstItTgtMat'] = False
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