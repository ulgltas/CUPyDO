#! /usr/bin/env python
# -*- coding: utf-8 -*-
# original name: fsi_StaticCylinder_Mtf_Pfem.py

import os
import sys

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


def getParameters(_p):
    p = {}
    p['nthreads'] = 1
    p['nDim'] = 2
    p['tollFSI'] = 1e-6
    p['dt'] = 0.1
    p['tTot'] = 1.0
    p['nFSIIterMax'] = 20
    p['timeIterTreshold'] = 0
    p['omegaMax'] = 0.5
    p['computationType'] = 'unsteady'
    p['saveFreqPFEM'] = 1
    p['mtfSaveAllFacs'] = False
    p.update(_p)
    return p


def main(_p, nogui):  # NB, the argument 'nogui' is specific to PFEM only!

    # --- Get FSI parameters ---#
    p = getParameters(_p)

    # --- Set up MPI --- #
    withMPI, comm, myid, numberPart = cupyutil.getMpi()
    rootProcess = 0
    #cupyutil.load(filePath, fileName, withMPI, comm, myid, numberPart)

    # --- Input parameters --- #
    cfd_file = 'cylinder_fluid'
    csd_file = 'cylinder_solid'

    # --- Initialize the fluid solver --- #
    import cupydoInterfaces.PfemInterface
    fluidSolver = cupydoInterfaces.PfemInterface.PfemSolver(
        cfd_file, 9, p['dt'])

    # --- This part is specific to PFEM ---
    fluidSolver.pfem.scheme.nthreads = p['nthreads']
    fluidSolver.pfem.scheme.savefreq = p['saveFreqPFEM']
    if nogui:
        fluidSolver.pfem.gui = None
    # ---

    cupyutil.mpiBarrier(comm)

    # --- Initialize the solid solver --- #
    solidSolver = None
    if myid == rootProcess:
        import cupydoInterfaces.MtfInterface
        solidSolver = cupydoInterfaces.MtfInterface.MtfSolver(
            csd_file, p['computationType'])

        # --- This part is specific to Metafor ---
        solidSolver.saveAllFacs = p['mtfSaveAllFacs']

    cupyutil.mpiBarrier(comm)

    # --- Initialize the FSI manager --- #
    manager = cupyman.Manager(fluidSolver, solidSolver,
                              p['nDim'], p['computationType'], comm)
    cupyutil.mpiBarrier()

    # --- Initialize the interpolator --- #
    interpolator = cupyinterp.MatchingMeshesInterpolator(
        manager, fluidSolver, solidSolver, comm)

    # --- Initialize the FSI criterion --- #
    criterion = cupycrit.DispNormCriterion(p['tollFSI'])
    cupyutil.mpiBarrier()

    # --- Initialize the FSI algorithm --- #
    algorithm = cupyalgo.AlgorithmBGSAitkenRelax(manager, fluidSolver,
                                                 solidSolver, interpolator,
                                                 criterion, p['nFSIIterMax'],
                                                 p['dt'], p['tTot'],
                                                 p['timeIterTreshold'], p['omegaMax'],
                                                 comm)

    # --- Launch the FSI computation --- #
    algorithm.run()

    # --- Check the results --- #

    # Check convergence and results
    if (algorithm.errValue > p['tollFSI']):
        print "\n\n" + "FSI residual = " + str(algorithm.errValue) + ", FSI tolerance = " + str(p['tollFSI'])
        raise Exception(ccolors.ANSI_RED +
                        "FSI algo failed to converge!" + ccolors.ANSI_RESET)

    # Read results from file
    with open("Node_6_POS.ascii", 'rb') as f:
        lines = f.readlines()
    result_1 = np.genfromtxt(lines[-1:], delimiter=None)
    with open("TOTAL_F_INT_Cylinder.ascii", 'rb') as f:
        lines = f.readlines()
    result_2 = np.genfromtxt(lines[-1:], delimiter=None)

    tests = CTests()
    tests.add(CTest('Mean nb of FSI iterations',
                    algorithm.getMeanNbOfFSIIt(), 3, 1, True))  # abs. tol. of 1
    tests.add(CTest('Y-coordinate Node 6',
                    result_1[1], 0.25, 1e-4, False))  # rel. tol. of 0.001%
    tests.add(CTest('Total force on cylinder',
                    result_2[1], -0.116, 1e-2, False))  # rel. tol. of 1%
    tests.run()

# -------------------------------------------------------------------
#    Run Main Program
# -------------------------------------------------------------------


# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':

    p = {}

    parser = OptionParser()
    parser.add_option("--nogui", action="store_true",
                      help="Specify if we need to use the GUI", dest="nogui", default=False)
    parser.add_option("-n", type="int", help="Number of threads",
                      dest="nthreads", default=1)

    (options, args) = parser.parse_args()

    nogui = options.nogui
    p['nthreads'] = options.nthreads

    main(p, nogui)
