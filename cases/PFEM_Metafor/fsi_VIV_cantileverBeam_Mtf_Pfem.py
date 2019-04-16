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
    p.update(_p)
    return p

def main(_p, nogui): # NB, the argument 'nogui' is specific to PFEM only!
    
    p = getParameters(_p)

    # --- Workspace set up --- #
    withMPI = False
    comm = None
    myid = 0
    numberPart = 0
    rootProcess = 0
    
    #cupyutil.load(filePath, fileName, withMPI, comm, myid, numberPart)
    
    cfd_file = 'VIV_cantileverBeam_air_Pfem'
    csd_file = 'VIV_cantileverBeam_beam_Mtf'
    nDim = 2
    tollFSI = 1e-6
    dt = 0.0005
    tTot = 10.
    nFSIIterMax = 10
    omegaMax = 1.0
    computationType = 'unsteady'
    
    # --- Initialize the fluid solver --- #
    import cupydoInterfaces.PfemInterface
    fluidSolver = cupydoInterfaces.PfemInterface.PfemSolver(cfd_file, 23, dt)
    
    # --- This part is specific to PFEM ---
    fluidSolver.pfem.scheme.nthreads = p['nthreads']
    if battery:
        fluidSolver.pfem.scheme.savefreq = int(tTot/dt)
    if nogui:
        fluidSolver.pfem.gui = None
    # ---
    
    cupyutil.mpiBarrier(comm)
    
    # --- Initialize the solid solver --- #
    solidSolver = None
    if myid == rootProcess:
        import cupydoInterfaces.MtfInterface
        solidSolver = cupydoInterfaces.MtfInterface.MtfSolver(csd_file, computationType)
        
        # --- This part is specific to Metafor ---
        if battery:
            solidSolver.saveAllFacs = False
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
    algorithm = cupyalgo.AlgorithmBGSAitkenRelax(manager, fluidSolver, solidSolver, interpolator, criterion, nFSIIterMax, dt, tTot, 0, omegaMax, comm)
    
    # --- Launch the FSI computation --- #
    algorithm.run()
    
# -------------------------------------------------------------------
#    Run Main Program
# -------------------------------------------------------------------

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':
    
    p = {}
    p['nthreads'] = nthreads
    
    parser=OptionParser()
    parser.add_option("--nogui", action="store_true",
                        help="Specify if we need to use the GUI", dest="nogui", default=False)
    parser.add_option("-n", type="int", help="Number of threads", dest="nthreads", default=1)


    (options, args)=parser.parse_args()
    
    nogui = options.nogui
    nthreads = options.nthreads
    
    main(p, nogui)