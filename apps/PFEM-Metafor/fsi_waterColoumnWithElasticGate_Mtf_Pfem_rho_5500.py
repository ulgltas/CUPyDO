import os, sys

filePath = os.path.abspath(os.path.dirname(sys.argv[0]))
fileName = os.path.splitext(os.path.basename(__file__))[0]
FSICouplerPath = 'D:/fsi'

sys.path.append(FSICouplerPath)

from math import *
from optparse import OptionParser
import FSICoupler

def main(battery, nogui): # NB, the two arguments 'battery' and 'nogui' are specific to PFEM problems involving FSI!
    
    # --- Workspace set up --- #
    withMPI = False
    comm = None
    myid = 0
    numberPart = 0
    rootProcess = 0
    
    FSICoupler.load(fileName, withMPI, comm, myid, numberPart)
    
    # --- Input parameters --- #
    cfd_file = 'waterColoumnWithElasticGate_water_Pfem'
    csd_file = 'waterColoumnWithElasticGate_gate_Mtf_rho_5500'
    nDim = 2
    tollFSI = 1e-6
    dt = 0.001
    tTot = 0.5
    nFSIIterMax = 10
    omegaMax = 0.5
    computationType = 'unsteady'
    
    # --- Initialize the fluid solver --- #
    import PfemInterface
    fluidSolver = PfemInterface.PfemSolver(cfd_file, 17, dt)
    
    # --- This part is specific to PFEM ---
    if battery:
        fluidSolver.pfem.scheme.savefreq = int(tTot/dt)
    if nogui:
        fluidSolver.pfem.gui = None
    # ---
    
    FSICoupler.mpiBarrier(comm)
    
    # --- Initialize the solid solver --- #
    solidSolver = None
    if myid == rootProcess:
        import MtfInterface
        solidSolver = MtfInterface.MtfSolver(csd_file, computationType)
        
        # --- This part is specific to Metafor ---
        if battery:
            solidSolver.saveAllFacs = False
        # ---
        
    FSICoupler.mpiBarrier(comm)
        
    # --- Initialize the FSI manager --- #
    manager = FSICoupler.Manager(fluidSolver, solidSolver, nDim, computationType, comm)
    FSICoupler.mpiBarrier()

    # --- Initialize the interpolator --- #
    interpolator = FSICoupler.MatchingMeshesInterpolator(manager, fluidSolver, solidSolver, comm)
    
    # --- Initialize the FSI criterion --- #
    criterion = FSICoupler.DispNormCriterion(tollFSI)
    FSICoupler.mpiBarrier()
    
    # --- Initialize the FSI algorithm --- #
    algorithm = FSICoupler.AlgortihmBGSAitkenRelax(manager, fluidSolver, solidSolver, interpolator, criterion, nFSIIterMax, dt, tTot, 0, omegaMax, comm)
    
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
    
    parser=OptionParser()
    parser.add_option("--nogui", action="store_true",
                        help="Specify if we need to use the GUI", dest="nogui", default=False)

    (options, args)=parser.parse_args()
    
    nogui = options.nogui
    
    if len(sys.argv) > 2:
        battery=sys.argv[1]
        main(battery, nogui)
    else:
        main(False, nogui)