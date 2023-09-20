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
    p['tollFSI'] = 1e-7
    p['dt'] = 1e-4
    p['tTot'] = 0.5
    p['nFSIIterMax'] = 20
    p['timeIterTreshold'] = 0
    p['omegaMax'] = 0.5
    p['aitkenCrit'] = 'min'
    p['regime'] = 'unsteady'
    p['testName'] = fileName
    p['cfdFile'] = 'waterColumnOnElasticColumn_water_Pfem'
    p['csdFile'] = 'waterColumnOnElasticColumn_elasticColumn_Mtf'
    p['betaFSI'] = 0.
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
    cfd_file = p['cfdFile']
    csd_file = p['csdFile']
    
    # --- Initialize the fluid solver --- #
    import cupydo.interfaces.Pfem as fItf
    fluidSolver = fItf.Pfem(cfd_file, 14, p['dt'])
    fluidSolver.pfem.pbl.betaFSI = p['betaFSI']
    
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
    algorithm = cupyalgo.AlgorithmBGSAitkenRelax(manager, fluidSolver, solidSolver, interpolator, criterion, p['nFSIIterMax'], p['dt'], p['tTot'], p['timeIterTreshold'], p['omegaMax'], comm)
    algorithm.atikenCrit = p['aitkenCrit']
    
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