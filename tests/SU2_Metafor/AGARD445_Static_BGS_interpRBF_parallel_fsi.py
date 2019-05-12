#! /usr/bin/env python
# -*- coding: utf8 -*-

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

import os, sys

filePath = os.path.abspath(os.path.dirname(__file__))
fileName = os.path.splitext(os.path.basename(__file__))[0]

from math import *
from optparse import OptionParser

import cupydo.utilities as cupyutil
import cupydo.manager as cupyman
import cupydo.interpolator as cupyinterp
import cupydo.criterion as cupycrit
import cupydo.algorithm as cupyalgo

import numpy as np
from cupydo.testing import *

def test(nogui, res, tol):

    # Read results from file
    with open("AerodynamicCoeff.ascii", 'rb') as f:
        lines = f.readlines()
    resultA = np.genfromtxt(lines[-1:], delimiter=None)
    with open("db_Field(TZ,RE)_GROUP_ID_180.ascii", 'rb') as f:
        lines = f.readlines()
    resultS1 = np.genfromtxt(lines[-1:], delimiter=None)
    with open("db_Field(TZ,RE)_GROUP_ID_181.ascii", 'rb') as f:
        lines = f.readlines()
    resultS2 = np.genfromtxt(lines[-1:], delimiter=None)

    # Check convergence and results
    # residual for test is 0.002
    #if (res > tol):
    #    print "\n\n" + "FSI residual = " + str(res) + ", FSI tolerance = " + str(tol)
    #    raise Exception(ccolors.ANSI_RED + "FSI algo failed to converge!" + ccolors.ANSI_RESET)
    tests = CTests()
    tests.add(CTest('Lift coefficient', resultA[2], 0.0537, 1e-1, False)) # rel. tol. of 10%
    tests.add(CTest('Drag coefficient', resultA[3], 0.00035, 1e-1, False)) # rel. tol. of 10%
    tests.add(CTest('Displacement (180, TZ)', resultS1[2], 0.0116, 1e-1, False)) # rel. tol. of 10%
    tests.add(CTest('Displacement (181, TZ)', resultS2[2], 0.0132, 1e-1, False)) # rel. tol. of 10%
    tests.run()

def getParameters(_p):
    # --- Input parameters --- #
    p = {}
    p['nDim'] = 3
    p['tollFSI'] = 1e-6
    p['dt'] = 0.0
    p['tTot'] = 0.05
    p['nFSIIterMax'] = 4
    p['rbfRadius'] = 0.3
    p['timeIterTreshold'] = -1
    p['omegaMax'] = 1.0
    p['computationType'] = 'steady'
    p['mtfSaveAllFacs'] = False
    p['nodalLoadsType'] = 'force'
    p['nZones_SU2'] = 0
    p.update(_p)
    return p 

def main(_p, nogui):

    # --- Get FSI parameters ---#
    p = getParameters(_p)

    # --- Set up MPI --- #
    withMPI, comm, myid, numberPart = cupyutil.getMpi()
    rootProcess = 0

    # --- Input parameters --- #
    cfd_file = '../../tests/SU2_Metafor/AGARD445_Static_SU2Conf.cfg'
    csd_file = 'AGARD445_Static_MetaforConf'

    # --- Initialize the fluid solver --- #
    import cupydo.interfaces.SU2 as fItf
    if comm != None:
        FluidSolver = fItf.SU2(cfd_file, p['nZones_SU2'], p['nDim'], p['computationType'], p['nodalLoadsType'], withMPI, comm)
    else:
        FluidSolver = fItf.SU2(cfd_file, p['nZones_SU2'], p['nDim'], p['computationType'], p['nodalLoadsType'], withMPI, 0)
    cupyutil.mpiBarrier(comm)

    # --- Initialize the solid solver --- #
    SolidSolver = None
    if myid == rootProcess:
        import cupydo.interfaces.Metafor as sItf
        SolidSolver = sItf.Metafor(csd_file, p['computationType'])
        SolidSolver.saveAllFacs = p['mtfSaveAllFacs']
    cupyutil.mpiBarrier(comm)

    # --- Initialize the FSI manager --- #
    manager = cupyman.Manager(FluidSolver, SolidSolver, p['nDim'], p['computationType'], comm)
    cupyutil.mpiBarrier()

    # --- Initialize the interpolator --- #
    interpolator = cupyinterp.RBFInterpolator(manager, FluidSolver, SolidSolver, p['rbfRadius'], comm)
    solverList = interpolator.getLinearSolvers()
    for ii in range(2):
        solverList[ii].setMaxNumberIterations(1000)
        solverList[ii].setPreconditioner("JACOBI")

    # --- Initialize the FSI criterion --- #
    criterion = cupycrit.DispNormCriterion(p['tollFSI'])
    cupyutil.mpiBarrier()

    # --- Initialize the FSI algorithm --- #
    algorithm = cupyalgo.AlgorithmBGSStaticRelax(manager, FluidSolver, SolidSolver, interpolator, criterion, p['nFSIIterMax'], p['dt'], p['tTot'], p['timeIterTreshold'], p['omegaMax'], comm)

    # --- Launch the FSI computation --- #
    algorithm.run()

    # --- Check the results --- #
    test(nogui, algorithm.errValue, p['tollFSI'])
  
    # --- Exit computation --- #
    del manager
    del criterion
    del FluidSolver
    del SolidSolver
    del interpolator
    del algorithm
    cupyutil.mpiBarrier(comm)
    return 0


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':

    p = {}
    
    parser=OptionParser()
    parser.add_option("--nogui", action="store_true",
                        help="Specify if we need to use the GUI", dest="nogui", default=False)
    parser.add_option("-n", type="int", help="Number of process", dest="nprocess", default=1) # not used

    (options, args)=parser.parse_args()
    
    nogui = options.nogui
    
    main(p, nogui)
