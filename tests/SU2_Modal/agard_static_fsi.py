#!/usr/bin/env python

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys

filePath = os.path.abspath(os.path.dirname(sys.argv[0]))
fileName = os.path.splitext(os.path.basename(__file__))[0]

from optparse import OptionParser	# use a parser for configuration
from math import *
import numpy as np

import cupydo.utilities as cupyutil
import cupydo.manager as cupyman
import cupydo.interpolator as cupyinterp
import cupydo.criterion as cupycrit
import cupydo.algorithm as cupyalgo

# -------------------------------------------------------------------
#    Script
# -------------------------------------------------------------------

def getParameters(_p):
    # --- Input parameters --- #
    p = {}
    p['nDim'] = 3
    p['tollFSI'] = 1e-6
    p['dt'] = 0.
    p['tTot'] = 0.
    p['nFSIIterMax'] = 5
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

    # --- Set up MPI and workspace --- #
    withMPI, comm, myid, numberPart = cupyutil.getMpi()
    rootProcess = 0
    cupyutil.load(fileName, withMPI, comm, myid, numberPart)

    # --- Input parameters --- #
    cfd_file = '../../tests/SU2_Modal/agard_static_fluid.cfg'
    csd_file = 'agard_solid'

    # --- Initialize the fluid solver --- #
    import cupydoInterfaces.SU2Interface
    if comm != None:
        fluidSolver = cupydoInterfaces.SU2Interface.SU2Solver(cfd_file, p['nZones_SU2'], p['nDim'], p['computationType'], p['nodalLoadsType'], withMPI, comm)
    else:
        fluidSolver = cupydoInterfaces.SU2Interface.SU2Solver(cfd_file, p['nZones_SU2'], p['nDim'], p['computationType'], p['nodalLoadsType'], withMPI, 0)
    cupyutil.mpiBarrier(comm)

    # --- Initialize modal interpreter --- #
    solidSolver = None
    if myid == rootProcess:
        import cupydoInterfaces.ModalInterface
        solidSolver = cupydoInterfaces.ModalInterface.ModalInterface(csd_file, p['computationType'])
    cupyutil.mpiBarrier(comm)

    # --- Initialize the FSI manager --- #
    manager = cupyman.Manager(fluidSolver, solidSolver, p['nDim'], p['computationType'], comm)
    cupyutil.mpiBarrier()

    # --- Initialize the interpolator --- #
    #interpolator = cupyinterp.MatchingMeshesInterpolator(manager, fluidSolver, solidSolver, comm)
    interpolator = cupyinterp.RBFInterpolator(manager, fluidSolver, solidSolver, p['rbfRadius'], comm)

    # --- Initialize the FSI criterion --- #
    criterion = cupycrit.DispNormCriterion(p['tollFSI'])
    cupyutil.mpiBarrier()

    # --- Initialize the FSI algorithm --- #
    algorithm = cupyalgo.AlgorithmBGSStaticRelax(manager, fluidSolver, solidSolver, interpolator, criterion, p['nFSIIterMax'], p['dt'], p['tTot'], p['timeIterTreshold'], p['omegaMax'], comm)
    
    # --- Launch the FSI computation --- #
    algorithm.run()

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
