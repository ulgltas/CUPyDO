import os, sys

filePath = os.path.abspath(os.path.dirname(sys.argv[0]))
fileName = os.path.splitext(os.path.basename(__file__))[0]
FSICouplerPath = os.path.join(os.path.dirname(__file__), '../../')

sys.path.append(FSICouplerPath)

from math import *
from optparse import OptionParser
import FSICoupler

def getParameters(_p):
    # --- Input parameters --- #
    p = {}
    p['nDim'] = 2
    p['tollFSI'] = 1e-6
    p['dt'] = 0.001
    p['tTot'] = 5.0
    p['nFSIIterMax'] = 20
    p['timeIterTreshold'] = 0
    p['omegaMax'] = 0.5
    p['computationType'] = 'unsteady'
    p['nbTimeToKeep'] = 0
    p['computeTangentMatrixBasedOnFirstIt'] = False
    p['saveFreqPFEM'] = 1
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
    
    FSICoupler.load(fileName, withMPI, comm, myid, numberPart)
    
    # --- Input parameters --- #
    cfd_file = 'waterColoumnFallWithFlexibleObstacle_water_Pfem'
    csd_file = 'waterColoumnFallWithFlexibleObstacle_obstacle_Mtf_E_1_0e6'
    
    # --- Initialize the fluid solver --- #
    import PfemInterface
    fluidSolver = PfemInterface.PfemSolver(cfd_file, 17, p['dt'])
    
    # --- This part is specific to PFEM ---
    fluidSolver.pfem.scheme.savefreq = p['saveFreqPFEM']
    if nogui:
        fluidSolver.pfem.gui = None
    # ---
    
    FSICoupler.mpiBarrier(comm)
    
    # --- Initialize the solid solver --- #
    solidSolver = None
    if myid == rootProcess:
        import MtfInterface
        solidSolver = MtfInterface.MtfSolver(csd_file, p['computationType'])
        
        # --- This part is specific to Metafor ---
        solidSolver.saveAllFacs = p['mtfSaveAllFacs']
        
    FSICoupler.mpiBarrier(comm)
        
    # --- Initialize the FSI manager --- #
    manager = FSICoupler.Manager(fluidSolver, solidSolver, p['nDim'], p['computationType'], comm)
    FSICoupler.mpiBarrier()

    # --- Initialize the interpolator --- #
    interpolator = FSICoupler.MatchingMeshesInterpolator(manager, fluidSolver, solidSolver, comm)
    
    # --- Initialize the FSI criterion --- #
    criterion = FSICoupler.DispNormCriterion(p['tollFSI'])
    FSICoupler.mpiBarrier()
    
    # --- Initialize the FSI algorithm --- #
    algorithm = FSICoupler.AlgortihmIQN_ILS(manager, fluidSolver, solidSolver, interpolator, criterion, p['nFSIIterMax'], p['dt'], p['tTot'], p['timeIterTreshold'], p['omegaMax'], p['nbTimeToKeep'], p['computeTangentMatrixBasedOnFirstIt'], comm)
    
    # --- Launch the FSI computation --- #
    algorithm.run()

    # --- Exit the fluid solver --- #
    fluidSolver.exit()

    # --- Exit the solid solver --- #
    if myid == rootProcess:
        solidSolver.exit()
    
    # --- Exit computation --- #
    FSICoupler.mpiBarrier(comm)
    return 0


# -------------------------------------------------------------------
#    Run Main Program
# -------------------------------------------------------------------

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':
    
    p = {}
    
    parser=OptionParser()
    parser.add_option("--nogui", action="store_true",
                        help="Specify if we need to use the GUI", dest="nogui", default=False)

    (options, args)=parser.parse_args()
    
    nogui = options.nogui
    
    main(p, nogui)