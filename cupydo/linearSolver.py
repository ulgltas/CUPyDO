#!/usr/bin/env python
# -*- coding: latin-1; -*-
#
# linearSolver.py
# Defines linear solver to be used with interface data.
# Authors : David THOMAS, Marco Lucio CERQUAGLIA, Romain BOMAN
#
# COPYRIGHT (C) University of Liège, 2017.

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import numpy as np
import scipy as sp
import scipy.sparse.linalg as splinalg

import ccupydo

np.set_printoptions(threshold=np.nan)

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
