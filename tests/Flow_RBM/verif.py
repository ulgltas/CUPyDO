#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from cupydo.testing import *

def staticAirfoil(nogui, res, tol):
    # Flow constant (defined in case_fluid)
    dynP = 0.5*100 # dynamic pressure
    alpha0 = np.radians(3) # initial angle of attack
    cRef = 1

    # RBM constant (defined in case_solid)
    k = 250 # vertical stiffness
    kappa = 25 # rotational stiffness

    # Airfoil slopes (measured from flow sovler)
    cl_alpha = 6.8
    cd_alpha = 0.085
    cm_alpha = 0.27 # must be measured from flexural axis posit. hard to calibrate but of crucial importance!

    # Pitch-plunge system of equations
    A = np.array([ [k, -dynP*cRef*(cl_alpha*np.cos(alpha0) - cd_alpha*np.sin(alpha0))],
                   [0, kappa - dynP*cRef*cRef*cm_alpha] ])
    b = np.array([0, kappa*alpha0])
    x = np.linalg.solve(A, b)

    # Display the solution
    if not nogui:
        print "Ref. lift coefficient: " + str(cl_alpha*x[1])
        print "Ref. vertical displacement: " + str(x[0])
        print "Ref. new angle of attack: " + str(np.degrees(x[1]))
        print "Ref. rotational displacement : " + str(np.degrees(x[1]-alpha0))

    # Read results from file
    resultS = np.genfromtxt("NativeHistory.dat", delimiter=None, skip_header=1)
    with open("FlowHistory.dat", 'rb') as f:
        lines = f.readlines()
    resultA = np.genfromtxt(lines[-1:], delimiter=None)

    # Check convergence and results
    if (res > tol):
        print "\n\n" + "FSI residual = " + str(res) + ", FSI tolerance = " + str(tol)
        raise Exception(ccolors.ANSI_RED + "FSI algo failed to converge!" + ccolors.ANSI_RESET)
    tests = CTests()
    tests.add(CTest('Lift coefficient', resultA[2], cl_alpha*x[1], 5e-2, False)) # rel. tol. of 5%
    tests.add(CTest('Vertical displacement', -resultS[2]/cRef, x[0]/cRef, 1e-2, True)) # abs. tol. of 1% of chord
    tests.add(CTest('Rotational displacement', np.degrees(resultS[3]), np.degrees(x[1]-alpha0), 5e-1, True)) # abs. tol. of .5°
    tests.run()
