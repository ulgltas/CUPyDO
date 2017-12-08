#!/usr/bin/env python
# -*- coding: latin-1; -*-
#
# FSICoupler.py
# Main file (Python core) of CUPyDO.
# Authors : David THOMAS, Marco Lucio CERQUAGLIA, Romain BOMAN
#
# COPYRIGHT (C) University of Liège, 2017.

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

from math import *
import numpy as np
import scipy as sp
from scipy import spatial
import scipy.sparse.linalg as splinalg
import os, os.path, sys, string
import time as tm

import traceback

import socket, fnmatch
import fsi_pyutils

import copy

import ccupydo

np.set_printoptions(threshold=np.nan)

# global vars (underscore prevent them to be imported with "from module import *")
_theModule  = None
_theWDir    = None # workspace directory
_theWDirRoot = os.getcwd()  # base directory du workspace

# ----------------------------------------------------------------------
#  Linear solver class
# ----------------------------------------------------------------------

class LinearSolver(ccupydo.CLinearSolver):
    """
    Define a parallel linear solver (serial compatible).
    Designed to be used with InterfaceData and InterfaceMatrix classes.
    Inherited public members :
        -solve()
    """

    def __init__(self, MatrixOperator, mpiComm=None):
        """
        Constructor.
        MatrixOperator is of type InterfaceMatrix
        """

        ccupydo.CLinearSolver.__init__(self, MatrixOperator)

        self.mpiComm = mpiComm

        if mpiComm == None:
            self.LinOperator = MatrixOperator.getMat()

    def solve(self, DataB, DataX):
        """
        Solve system MatrixOperator*VecX = VecB.
        VecX and VecB are InterfaceData types.
        """

        if self.mpiComm != None:
            ccupydo.CLinearSolver.solve(self, DataB, DataX)
        else:
            dim = DataB.getDim()
            for iDim in range(dim):
                XX = splinalg.spsolve(self.LinOperator, DataB.getData(iDim))
                DataX.setData(iDim, XX)
