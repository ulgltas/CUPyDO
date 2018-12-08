#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from cupydo.testing import *

def staticAirfoil(nogui):
    # Flow constant (defined in case_fluid)
    dynP = 0.5*100 # dynamic pressure
    alpha0 = np.radians(3) # initial angle of attack
    cRef = 1

    # RBM constant (defined in case_solid)
    k = 500 # vertical stiffness
    kappa = 50 # rotational stiffness

    # Airfoil slopes (measured from flow sovler)
    cl_alpha = 6.8
    cd_alpha = 0.085
    cm_alpha = 0.25 # must be measured from flexural axis posit. hard to calibrate but of crucial importance!

    # Pitch-plunge system of equations
    A = np.array([ [k, -dynP*cRef*(cl_alpha*np.cos(alpha0) - cd_alpha*np.sin(alpha0))],
                   [0, kappa - dynP*cRef*cRef*cm_alpha] ])
    b = np.array([0, kappa*alpha0])
    x = np.linalg.solve(A, b)

    # Display the solution
    if not nogui:
        print "Ref. vertical displacement: " + str(x[0])
        print "Ref. new angle of attack: " + str(np.degrees(x[1]))
        print "Ref. rotational displacement : " + str(np.degrees(x[1]-alpha0))

    # Check the solution
    res = np.genfromtxt("NativeHistory.dat", delimiter=None, skip_header=1) # import results by reading file
    tests = CTests()
    tests.add(CTest('Vertical displacement', -res[2]/cRef, x[0]/cRef, 5e-3, True)) # abs. tol. of .5% of chord
    tests.add(CTest('Rotational displacement', np.degrees(res[3]), np.degrees(x[1]-alpha0), 1e-1, True)) # abs. tol. of .1°
    tests.run()
