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

def getParameters(_p):
    # --- Input parameters --- #
    p = {}
    p['nDim'] = 2
    p['tollFSI'] = 1e-6
    p['dt'] = 0.0005
    p['tTot'] = 10.0
    p['nFSIIterMax'] = 20
    p['timeIterTreshold'] = 0
    p['omegaMax'] = 0.5
    p['regime'] = 'unsteady'
    p['nbTimeToKeep'] = 30
    p['computeTangentMatrixBasedOnFirstIt'] = False
    p['QR_filter'] = 'Haelterman'
    p['tollQR'] = 1.0e-1
    p['extractor'] = None
    p.update(_p)
    return p

def main(_p):
    
    p = getParameters(_p)

    # --- Workspace set up --- #
    withMPI = False
    comm = None
    myid = 0
    numberPart = 0
    rootProcess = 0
    
    
    
    # --- Input parameters --- #
    cfd_file = 'elasticContainerFilling_fluid_mu_50_Pfem'
    csd_file = 'elasticContainerFilling_container_Mtf'
    
    # --- Initialize the fluid solver --- #
    import cupydo.interfaces.Pfem as fItf
    fluidSolver = fItf.Pfem(cfd_file, 5, p['dt'])
    
    cupyutil.mpiBarrier(comm)
    
    # --- Initialize the solid solver --- #
    solidSolver = None
    if myid == rootProcess:
        import cupydo.interfaces.Metafor as sItf
        solidSolver = sItf.Metafor(csd_file, p['regime'])
        
    cupyutil.mpiBarrier(comm)
        
    # --- Initialize the FSI manager --- #
    manager = cupyman.Manager(fluidSolver, solidSolver, p['nDim'], p['regime'], comm)
    cupyutil.mpiBarrier()

    # --- Initialize the interpolator --- #
    interpolator = cupyinterp.MatchingMeshesInterpolator(manager, fluidSolver, solidSolver, comm)
    
    # --- Initialize the FSI criterion --- #
    criterion = cupycrit.DispNormCriterion(p['tollFSI'])
    cupyutil.mpiBarrier()
    
    # --- Initialize the FSI algorithm --- #
    algorithm = cupyalgo.AlgorithmIQN_ILS(manager, fluidSolver, solidSolver, interpolator, criterion, p['nFSIIterMax'], p['dt'], p['tTot'], p['timeIterTreshold'], p['omegaMax'], p['nbTimeToKeep'], p['computeTangentMatrixBasedOnFirstIt'], comm)
    algorithm.qrFilter = p['QR_filter']
    algorithm.tollQR = p['tollQR']
    
    # --- Launch the FSI computation --- #
    algorithm.run()

# -------------------------------------------------------------------
#    Run Main Program
# -------------------------------------------------------------------

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':
    
    p = {}
    parser=OptionParser()
    parser.add_option("-k", type="int", help="Number of threads", dest="nthreads", default=1)
    (options, args)=parser.parse_args()
    p['nthreads'] = options.nthreads
    main(p)