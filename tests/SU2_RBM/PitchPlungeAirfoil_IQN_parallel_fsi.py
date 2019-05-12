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

    # Check convergence and results
    if (res > tol):
        print "\n\n" + "FSI residual = " + str(res) + ", FSI tolerance = " + str(tol)
        raise Exception(ccolors.ANSI_RED + "FSI algo failed to converge!" + ccolors.ANSI_RESET)
    tests = CTests()
    tests.add(CTest('Lift coefficient', resultA[2], 0.245, 1e-1, False)) # rel. tol. of 10%
    tests.add(CTest('Drag coefficient', resultA[3], 0.0015, 1e-1, False)) # rel. tol. of 10%
    tests.run()

def getParameters(_p):
    # --- Input parameters --- #
    p = {}
    p['nDim'] = 2
    p['tollFSI'] = 1e-8
    p['dt'] = 0.001
    p['tTot'] = 0.005
    p['nFSIIterMax'] = 6
    p['timeIterTreshold'] = 0
    p['omegaMax'] = 1.0
    p['nbTimeToKeep'] = 0
    p['computeTangentMatrixBasedOnFirstIt'] = False
    p['computationType'] = 'unsteady'
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

    cfd_file = '../../tests/SU2_RBM/PitchPlungeAirfoil_BGS_parallel_SU2Conf.cfg'
    csd_file = '../../tests/SU2_RBM/PitchPlungeAirfoil_BGS_parallel_RBMConf.cfg'

    # --- Initialize the fluid solver --- #
    import cupydo.interfaces.SU2 as fItf
    if comm != None:
      fluidSolver = fItf.SU2(cfd_file, p['nZones_SU2'], p['nDim'], p['computationType'], p['nodalLoadsType'], withMPI, comm)
    else:
      fluidSolver = fItf.SU2(cfd_file, p['nZones_SU2'], p['nDim'], p['computationType'], p['nodalLoadsType'], withMPI, 0)
  
    cupyutil.mpiBarrier(comm)

    # --- Initialize the solid solver --- #
    solidSolver = None
    if myid == rootProcess:
      import cupydo.interfaces.RBMI as sItf
      solidSolver = sItf.RBMI(csd_file, p['computationType'])

    cupyutil.mpiBarrier(comm)

    # --- Initialize the FSI manager --- #
    manager = cupyman.Manager(fluidSolver, solidSolver, p['nDim'], p['computationType'], comm)
    cupyutil.mpiBarrier()

    # --- Initialize the interpolator --- #
    interpolator = cupyinterp.MatchingMeshesInterpolator(manager, fluidSolver, solidSolver, comm)

    # --- Initialize the FSI criterion --- #
    criterion = cupycrit.DispNormCriterion(p['tollFSI'])
    cupyutil.mpiBarrier()

    # --- Initialize the FSI algorithm --- #
    algorithm = cupyalgo.AlgorithmIQN_ILS(manager, fluidSolver, solidSolver, interpolator, criterion, p['nFSIIterMax'], p['dt'], p['tTot'], p['timeIterTreshold'], p['omegaMax'], p['nbTimeToKeep'], p['computeTangentMatrixBasedOnFirstIt'], comm)

    # --- Launch the FSI computation --- #
    algorithm.run()

    # --- Check the results --- #
    test(nogui, algorithm.errValue, p['tollFSI'])
  
    # --- Exit computation --- #
    del manager
    del criterion
    del fluidSolver
    del solidSolver
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
