import os, sys

filePath = os.path.abspath(os.path.dirname(sys.argv[0]))
fileName = os.path.splitext(os.path.basename(__file__))[0]
FSICouplerPath = os.path.join(os.path.dirname(__file__), '../../')

sys.path.append(FSICouplerPath)

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
    p['nDim'] = 2
    p['tollFSI'] = 1e-7
    p['dt'] = 1e-5
    p['tTot'] = 0.5
    p['nFSIIterMax'] = 20
    p['timeIterTreshold'] = 0
    p['omegaMax'] = 1.0
    p['computationType'] = 'unsteady'
    p['saveFreqPFEM'] = 100
    p['mtfSaveAllFacs'] = True
    p['testName'] = fileName
    p['cfdFile'] = 'waterColumnOnElasticColumn_water_Pfem'
    p['csdFile'] = 'waterColumnOnElasticColumn_elasticColumn_Mtf_rho_1500_E_2_3e5_nu_0_4'
    p['betaFSI'] = 0.
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
    
    cupyutil.load(p['testName'], withMPI, comm, myid, numberPart)
    
    # --- Input parameters --- #
    cfd_file = p['cfdFile']
    csd_file = p['csdFile']
    
    # --- Initialize the fluid solver --- #
    import cupydoInterfaces.PfemInterface
    fluidSolver = cupydoInterfaces.PfemInterface.PfemSolver(cfd_file, 14, p['dt'])
    fluidSolver.pfem.pbl.betaFSI = p['betaFSI']
    
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
        solidSolver = cupydoInterfaces.MtfInterface.MtfSolver(csd_file, p['computationType'])
        
        # --- This part is specific to Metafor ---
        solidSolver.saveAllFacs = p['mtfSaveAllFacs']
        
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
    algorithm = cupyalgo.AlgorithmBGSStaticRelax(manager, fluidSolver, solidSolver, interpolator, criterion, p['nFSIIterMax'], p['dt'], p['tTot'], p['timeIterTreshold'], p['omegaMax'], comm)
    
    # --- Launch the FSI computation --- #
    algorithm.run()

# -------------------------------------------------------------------
#    Run Main Program
# -------------------------------------------------------------------

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':
    
    p = {}
    
    parser=OptionParser()
    parser.add_option("--nogui", action="store_true",
                        help="Specify if we need to use the GUI", dest="nogui", default=False)
    parser.add_option("--nthreads", type="int", help="Number of threads", dest="nthreads", default=1)


    (options, args)=parser.parse_args()
    
    nogui = options.nogui
    nthreads = options.nthreads
    
    p['nthreads'] = nthreads
    
    main(p, nogui)