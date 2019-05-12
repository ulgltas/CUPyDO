#!/usr/bin/env python
# -*- coding: utf8 -*-

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys

filePath = os.path.abspath(os.path.dirname(__file__))
fileName = os.path.splitext(os.path.basename(__file__))[0]

from optparse import OptionParser	# use a parser for configuration
from math import *
import numpy as np

import cupydo.utilities as cupyutil
import cupydo.manager as cupyman
import cupydo.interpolator as cupyinterp
import cupydo.criterion as cupycrit
import cupydo.algorithm as cupyalgo
from cupydo.testing import *

# -------------------------------------------------------------------
#    Script
# -------------------------------------------------------------------

def test(res, tol):

    # Read results from file
    with open("AerodynamicCoeff.ascii", 'rb') as f:
        lines = f.readlines()
    resultA = np.genfromtxt(lines[-1:], delimiter=None)
    with open("NodalDisplacement.dat", 'rb') as f:
        lines = f.readlines()
    resultS = np.genfromtxt(lines[-1:], delimiter=None)

    # Check convergence and results
    if (res > tol):
        print "\n\n" + "FSI residual = " + str(res) + ", FSI tolerance = " + str(tol)
        raise Exception(ccolors.ANSI_RED + "FSI algo failed to converge!" + ccolors.ANSI_RESET)
    tests = CTests()
    tests.add(CTest('Lift coefficient', resultA[2], 0.0537, 1e-1, False)) # rel. tol. of 10%
    tests.add(CTest('Drag coefficient', resultA[3], 0.00035, 1e-1, False)) # rel. tol. of 10%
    tests.add(CTest('LE Displacement (16, z)', resultS[4], 0.0116, 1e-1, False)) # rel. tol. of 10%
    tests.add(CTest('TE Displacement (13808, z)', resultS[7], 0.0132, 1e-1, False)) # rel. tol. of 10%
    tests.run()

def getParameters(_p):
    # --- Input parameters --- #
    p = {}
    p['nDim'] = 3
    p['tollFSI'] = 5e-3
    p['dt'] = 0.
    p['tTot'] = 0.
    p['nFSIIterMax'] = 6
    p['rbfRadius'] = 1.0
    p['timeIterTreshold'] = -1
    p['omegaMax'] = 1.0
    p['computationType'] = 'steady'
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
    cfd_file = os.path.join(filePath,'staticAgard_fluid.cfg')
    csd_file = 'agard_solid'

    # --- Initialize the fluid solver --- #
    import cupydo.interfaces.SU2 as fItf
    if comm != None:
        fluidSolver = fItf.SU2(cfd_file, p['nZones_SU2'], p['nDim'], p['computationType'], p['nodalLoadsType'], withMPI, comm)
    else:
        fluidSolver = fItf.SU2(cfd_file, p['nZones_SU2'], p['nDim'], p['computationType'], p['nodalLoadsType'], withMPI, 0)
    cupyutil.mpiBarrier(comm)

    # --- Initialize modal solver --- #
    solidSolver = None
    if myid == rootProcess:
        import cupydo.interfaces.Modal as sItf
        solidSolver = sItf.Modal(csd_file, p['computationType'])
    cupyutil.mpiBarrier(comm)

    # --- Initialize the FSI manager --- #
    manager = cupyman.Manager(fluidSolver, solidSolver, p['nDim'], p['computationType'], comm)
    cupyutil.mpiBarrier()

    # --- Initialize the interpolator --- #
    interpolator = cupyinterp.RBFInterpolator(manager, fluidSolver, solidSolver, p['rbfRadius'], comm)

    # --- Initialize the FSI criterion --- #
    criterion = cupycrit.DispNormCriterion(p['tollFSI'])
    cupyutil.mpiBarrier()

    # --- Initialize the FSI algorithm --- #
    algorithm = cupyalgo.AlgorithmBGSStaticRelax(manager, fluidSolver, solidSolver, interpolator, criterion, p['nFSIIterMax'], p['dt'], p['tTot'], p['timeIterTreshold'], p['omegaMax'], comm)
    
    # --- Launch the FSI computation --- #
    algorithm.run()

    # --- Check the results --- #
    test(algorithm.errValue, p['tollFSI'])

    # --- Exit computation --- #
    del manager
    del criterion
    del fluidSolver
    del solidSolver
    del interpolator
    del algorithm
    cupyutil.mpiBarrier(comm)

    # eof
    print ''

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    p = {}
    
    parser=OptionParser()
    parser.add_option("--nogui", action="store_true",
                        help="Specify if we need to use the GUI", dest="nogui", default=False)
    parser.add_option("-n", type="int", help="Number of process", dest="nprocess", default=1) # not used

    (options, args)=parser.parse_args()
    
    nogui = options.nogui
    
    main(p, nogui)
