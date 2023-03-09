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
    p['nthreads'] = 1
    p['mtfSaveAllFacs'] = True
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
    
    
    
    cfd_file = 'VIV_cantileverBeam_air_Pfem'
    csd_file = 'VIV_cantileverBeam_beam_Mtf_EAS'
    nDim = 2
    tollFSI = 1e-6
    dt = 0.0005
    tTot = 10.
    nFSIIterMax = 10
    omegaMax = 1.0
    computationType = 'unsteady'
    
    # --- Initialize the fluid solver --- #
    import cupydo.interfaces.Pfem as fItf
    fluidSolver = fItf.Pfem(cfd_file, 23, dt)
    
    # --- This part is specific to PFEM ---
    fluidSolver.pfem.scheme.nthreads = p['nthreads']
    fluidSolver.pfem.scheme.savefreq = p['saveFreqPFEM']
    # ---
    
    cupyutil.mpiBarrier(comm)
    
    # --- Initialize the solid solver --- #
    solidSolver = None
    if myid == rootProcess:
        import cupydo.interfaces.Metafor as sItf
        solidSolver = sItf.Metafor(csd_file, computationType)
    # ---
        
    cupyutil.mpiBarrier(comm)
        
    # --- Initialize the FSI manager --- #
    manager = cupyman.Manager(fluidSolver, solidSolver, nDim, computationType, comm)
    cupyutil.mpiBarrier()

    # --- Initialize the interpolator --- #
    interpolator = cupyinterp.MatchingMeshesInterpolator(manager, fluidSolver, solidSolver, comm)
    
    # --- Initialize the FSI criterion --- #
    criterion = cupycrit.DispNormCriterion(tollFSI)
    cupyutil.mpiBarrier()
    
    # --- Initialize the FSI algorithm --- #
    algorithm = cupyalgo.AlgorithmBGSStaticRelax(manager, fluidSolver, solidSolver, interpolator, criterion, nFSIIterMax, dt, tTot, 0, omegaMax, comm)
    
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