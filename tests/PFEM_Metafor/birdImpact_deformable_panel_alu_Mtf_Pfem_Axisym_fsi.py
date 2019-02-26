import os, sys

filePath = os.path.abspath(os.path.dirname(sys.argv[0]))
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
    # --- Input parameters --- #
    p = {}
    p['nthreads'] = 1
    p['U0'] = 100
    p['N'] = 10
    p['R'] = 0.01
    p['d'] = 0.0025 # 2.5*(R/N)
    p['nDim'] = 2
    p['tollFSI'] = 1e-6
    p['dt'] = 1.5e-6
    p['tTot'] = 1e-4 # 40*((4*R)/U0 + d/U0)
    p['nFSIIterMax'] = 20
    p['timeIterTreshold'] = 0
    p['omegaMax'] = 0.5
    p['computationType'] = 'unsteady'
    p['saveFreqPFEM'] = 10
    p['mtfSaveAllFacs'] = False
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
    
    cupyutil.load(fileName, withMPI, comm, myid, numberPart)
    
    cfd_file = 'birdImpact_deformable_panel_bird_Pfem_Axisym'
    csd_file = 'birdImpact_deformable_panel_panel_alu_Mtf_Axisym'
    
    # --- Initialize the fluid solver --- #
    import cupydoInterfaces.PfemInterface
    fluidSolver = cupydoInterfaces.PfemInterface.PfemSolver(cfd_file, 13, p['dt'])
    
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
    algorithm = cupyalgo.AlgorithmBGSAitkenRelax(manager, fluidSolver, solidSolver, interpolator, criterion, p['nFSIIterMax'], p['dt'], p['tTot'], p['timeIterTreshold'], p['omegaMax'], comm)
    
    # --- Launch the FSI computation --- #
    algorithm.run()
    
    # --- Check the results --- #
    
    # Check convergence and results
    if (algorithm.errValue > p['tollFSI']):
        print "\n\n" + "FSI residual = " + str(algorithm.errValue) + ", FSI tolerance = " + str(p['tollFSI'])
        raise Exception(ccolors.ANSI_RED + "FSI algo failed to converge!" + ccolors.ANSI_RESET)
    
    # Read results from file
    with open("db_Field(TY,RE)_GROUP_ID_18.ascii", 'rb') as f:
        lines = f.readlines()
    result_1 = np.genfromtxt(lines[-1:], delimiter=None)
    
    tests = CTests()
    tests.add(CTest('Mean nb of FSI iterations', algorithm.getMeanNbOfFSIIt(), 4, 1, True)) # abs. tol. of 1
    tests.add(CTest('Y-displacement panel center', result_1[2], -0.002028, 1e-2, False)) # rel. tol. of 1%
    tests.run()

# -------------------------------------------------------------------
#    Run Main Program
# -------------------------------------------------------------------

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':
    
    p = {}
    
    parser=OptionParser()
    parser.add_option("--nogui", action="store_true",
                        help="Specify if we need to use the GUI", dest="nogui", default=False)
    parser.add_option("-n", type="int", help="Number of threads", dest="nthreads", default=1)
    
    
    (options, args)=parser.parse_args()
    
    nogui = options.nogui
    p['nthreads'] = options.nthreads
    
    main(p, nogui)