#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from cupydo.testing import *

def staticCylinder(nogui, res, tol):

	# Check convergence and results
    if (res > tol):
        print "\n\n" + "FSI residual = " + str(res) + ", FSI tolerance = " + str(tol)
        raise Exception(ccolors.ANSI_RED + "FSI algo failed to converge!" + ccolors.ANSI_RESET)

def birdImpact_steel(nogui, res, tol):
    
    # Read results from file
    with open("tf_THERMODYN_TRAV_FINT_GROUP_ID_17.ascii", 'rb') as f:
        lines = f.readlines()
    resultA = np.genfromtxt(lines[-1:], delimiter=None)
    with open("db_Field(TY,RE)_GROUP_ID_17.ascii", 'rb') as f:
        lines = f.readlines()
    resultS = np.genfromtxt(lines[-1:], delimiter=None)

    # Check convergence and results
    if (res > tol):
        print "\n\n" + "FSI residual = " + str(res) + ", FSI tolerance = " + str(tol)
        raise Exception(ccolors.ANSI_RED + "FSI algo failed to converge!" + ccolors.ANSI_RESET)
    tests = CTests()
    tests.add(CTest('Internal ? (17, FINT)', resultA[2], 4.9, 1e-1, False)) # rel. tol. of 10%
    tests.add(CTest('Displacement (17, TY)', resultS[2], -0.0007, 1e-1, False)) # rel. tol. of 10%
    tests.run()
