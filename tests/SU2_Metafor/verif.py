#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from cupydo.testing import *

def CantileverSquareChannel(nogui, res, tol):
    
    # Read results from file
    with open("AerodynamicCoeff.ascii", 'rb') as f:
        lines = f.readlines()
    resultA = np.genfromtxt(lines[-1:], delimiter=None)

    # Check convergence and results
    if (res > tol):
        print "\n\n" + "FSI residual = " + str(res) + ", FSI tolerance = " + str(tol)
        raise Exception(ccolors.ANSI_RED + "FSI algo failed to converge!" + ccolors.ANSI_RESET)
    tests = CTests()
    tests.add(CTest('Lift coefficient', resultA[2], 0.00095, 1e-1, False)) # rel. tol. of 10%
    tests.add(CTest('Drag coefficient', resultA[3], 2.65, 1e-1, False)) # rel. tol. of 10%
    tests.add(CTest('Displacement (104, TY)', resultS[2], 2.5e-8, 1e-1, False)) # rel. tol. of 10%
    tests.run()

def AGARD445_Static(nogui, res, tol):

        # Read results from file
    with open("AerodynamicCoeff.ascii", 'rb') as f:
        lines = f.readlines()
    resultA = np.genfromtxt(lines[-1:], delimiter=None)

    # Check convergence and results
    if (res > tol):
        print "\n\n" + "FSI residual = " + str(res) + ", FSI tolerance = " + str(tol)
        raise Exception(ccolors.ANSI_RED + "FSI algo failed to converge!" + ccolors.ANSI_RESET)
    tests = CTests()
    tests.add(CTest('Lift coefficient', resultA[2], 0.053, 1e-1, False)) # rel. tol. of 10%
    tests.add(CTest('Drag coefficient', resultA[3], 0.00035, 1e-1, False)) # rel. tol. of 10%
    tests.add(CTest('Displacement (180, TZ)', resultS[2], 0.011, 1e-1, False)) # rel. tol. of 10%
    tests.add(CTest('Displacement (181, TZ)', resultS[3], 0.013, 1e-1, False)) # rel. tol. of 10%
    tests.run()
