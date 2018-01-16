''' 

Copyright 2018 University of Liège

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 

'''

#!/usr/bin/env python
# -*- coding: latin-1; -*-
#
# linearSolver.py
# Defines linear solver to be used with interface data.
# Authors : David THOMAS, Marco Lucio CERQUAGLIA, Romain BOMAN

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
