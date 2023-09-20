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
    p['dt'] = 0.001
    p['tTot'] = 0.5
    p['nFSIIterMax'] = 20
    p['timeIterTreshold'] = 0
    p['omegaMax'] = 0.5
    p['regime'] = 'unsteady'
    p['nbTimeToKeep'] = 0
    p['computeTangentMatrixBasedOnFirstIt'] = False
    p['tollQR'] = 1e-3
    p['testName'] = fileName
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
    
    cupyutil.load(p['testName'], withMPI, comm, myid, numberPart)
    
    # --- Input parameters --- #
    cfd_file = 'waterColoumnWithElasticGate_water_Pfem_fine'
    csd_file = 'waterColoumnWithElasticGate_gate_Mtf_rho_1100_fine_gravity'
    
    # --- Initialize the fluid solver --- #
    import cupydo.interfaces.Pfem as fItf
    fluidSolver = fItf.Pfem(cfd_file, 17, p['dt'])
    
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