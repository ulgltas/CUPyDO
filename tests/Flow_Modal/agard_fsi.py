#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys

filePath = os.path.abspath(os.path.dirname(sys.argv[0]))
fileName = os.path.splitext(os.path.basename(__file__))[0]

from math import *
from optparse import OptionParser

import cupydo.utilities as cupyutil
import cupydo.manager as cupyman
import cupydo.interpolator as cupyinterp
import cupydo.criterion as cupycrit
import cupydo.algorithm as cupyalgo

from cupydo.testing import *
import numpy as np

def test(res, tol):
    # Read results from file
    with open("FlowHistory.dat", 'rb') as f:
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
    tests.add(CTest('Lift coefficient', resultA[2], 0.0537, 1e-2, True)) # abs. tol
    tests.add(CTest('Drag coefficient', resultA[3], 0.00045, 1e-4, True))
    tests.add(CTest('LE Displacement (16, z)', resultS[4], 0.0116, 1e-1, False)) # rel. tol. of 10%
    tests.add(CTest('TE Displacement (13808, z)', resultS[7], 0.0132, 1e-1, False)) # rel. tol. of 10%
    tests.run()

def getParameters(_p):
    # --- Input parameters --- #
    p = {}
    p['nthreads'] = 1
    p['nDim'] = 3
    p['tollFSI'] = 1e-5
    p['dt'] = 0.
    p['tTot'] = 0.
    p['nFSIIterMax'] = 50
    p['timeIterTreshold'] = -1
    p['RBFradius'] = 1.
    p['omegaMax'] = 1.0
    p['nbTimeToKeep'] = 0
    p['computeTangentMatrixBasedOnFirstIt'] = False
    p['computationType'] = 'steady'
    p.update(_p)
    return p

def main(_p, nogui):
    
    # --- Get FSI parameters ---#
    p = getParameters(_p)

    # --- Set up MPI and workspace --- #
    withMPI, comm, myid, numberPart = cupyutil.getMpi()
    rootProcess = 0
    cupyutil.load(fileName, withMPI, comm, myid, numberPart)
    
    # --- Input files --- #
    cfd_module = fileName[:-3] + "fluid"
    csd_module = fileName[:-3] + "solid"
    
    # --- Initialize the fluid solver --- #
    import cupydoInterfaces.FlowInterface
    fluidSolver = cupydoInterfaces.FlowInterface.Flow(cfd_module, p['nthreads'])
    
    cupyutil.mpiBarrier(comm)
    
    # --- Initialize modal solver --- #
    solidSolver = None
    if myid == rootProcess:
        import cupydoInterfaces.ModalInterface
        solidSolver = cupydoInterfaces.ModalInterface.ModalInterface(csd_module, p['computationType'])
    cupyutil.mpiBarrier(comm)
        
    # --- Initialize the FSI manager --- #
    manager = cupyman.Manager(fluidSolver, solidSolver, p['nDim'], p['computationType'], comm)
    cupyutil.mpiBarrier()

    # --- Initialize the interpolator --- #
    interpolator = cupyinterp.RBFInterpolator(manager, fluidSolver, solidSolver, p['RBFradius'], comm)
    
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
#    Run Main Program
# -------------------------------------------------------------------

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':
    
    p = {}
    
    parser=OptionParser()
    parser.add_option("--nogui", action="store_true",
                        help="Specify if we need to use the GUI", dest="nogui", default=False)
    parser.add_option("-n", type="int", help="Number of threads", dest="nthreads", default=1)
    
    
    (options, args)=parser.parse_args()
    
    nogui = options.nogui
    p['nthreads'] = options.nthreads
    
    main(p, nogui)