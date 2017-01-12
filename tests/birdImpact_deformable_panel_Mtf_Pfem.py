import os, sys

filePath = os.path.abspath(os.path.dirname(sys.argv[0]))
fileName = os.path.splitext(os.path.basename(__file__))[0]
fsiPath = '..'
sys.path.append(fsiPath)

import time as timer
from math import *
from optparse import OptionParser
import fsi
import PfemInterface
import MtfInterface

def main():
    
  fsi.load(fileName)
    
# --- Get the FSI conig file name form the command line options --- #
  parser=OptionParser()
  parser.add_option("-f", "--file",       dest="filename",
                      help="read config from FILE", metavar="FILE")
  parser.add_option("--parallel", action="store_true",
                      help="Specify if we need to initialize MPI", dest="with_MPI", default=False)

  (options, args)=parser.parse_args()

  if options.with_MPI == True:
    from mpi4py import MPI  # MPI is initialized from now by python and can be continued in C++ !
    comm = MPI.COMM_WORLD
    myid = comm.Get_rank()
    numberPart = comm.Get_size()
    have_MPI = True
  else:
    comm = None
    myid = 0
    numberPart = 1
    have_MPI = False

  rootProcess = 0

  # --- Input parameters --- #
  
  U0 = 100
  N = 10
  R = 0.01
  d = 2.5*(R/N)
    
  CFD_file = 'birdImpact_deformable_panel_bird_Pfem'
  CSD_file = 'birdImpact_deformable_panel_panel_steel_Mtf'
  nDim = 2
  FSITolerance = 1e-6
  deltaT = 2e-6
  totTime = 100*deltaT # 40*((4*R)/U0 + d/U0)
  nbFSIIterMax = 20
  omegaMax = 0.5
  computationType = 'unsteady'

  # --- Start timer --- #
  startTime = timer.time()

  # --- Initialize the fluid solver --- #
  fsi.MPIPrint('\n***************************** Initializing fluid solver *****************************', comm)
  # import PfemInterface
  FluidSolver = PfemInterface.PfemSolver(CFD_file, 13, deltaT)
  fsi.MPIBarrier(comm)

  # --- Initialize the solid solver --- #
  fsi.MPIPrint('\n***************************** Initializing solid solver *****************************', comm)
  SolidSolver = None
  if myid == rootProcess:
    # import MtfInterface
    SolidSolver = MtfInterface.MtfSolver(CSD_file)
  fsi.MPIBarrier(comm)

  # --- Initialize the FSI manager --- #
  fsi.MPIPrint('\n***************************** Initializing FSI interface *****************************', comm)
  manager = fsi.Manager(FluidSolver, SolidSolver, nDim, computationType, comm)
  fsi.MPIBarrier()

  # --- Initialize the interpolator --- #
  fsi.MPIPrint('\n***************************** Initializing FSI interpolator *****************************', comm)
  interpolator = fsi.MatchingMeshesInterpolator(manager, FluidSolver, SolidSolver, comm)

  # --- Initialize the FSI criterion --- #
  criterion = fsi.DispNormCriterion(FSITolerance)
  fsi.MPIBarrier()

  # --- Initialize the FSI algorithm --- #
  fsi.MPIPrint('\n***************************** Initializing FSI algorithm *****************************', comm)
  algorithm = fsi.AlgortihmBGSAitkenRelax(manager, FluidSolver, SolidSolver, interpolator, criterion, nbFSIIterMax, deltaT, totTime, 0, omegaMax, comm)

  # --- Launch the FSI computation --- #
  algorithm.run()

  # --- Exit the fluid solver --- #
  FluidSolver.exit()

  # --- Exit the solid solver --- #
  if myid == rootProcess:
    SolidSolver.exit()

  # --- Stop timer --- #
  stopTime = timer.time()
  elapsedTime = stopTime-startTime
  
  # --- Exit computation --- #
  del manager
  del criterion
  del FluidSolver
  del SolidSolver
  del interpolator
  del algorithm
  fsi.MPIBarrier(comm)
  fsi.MPIPrint("\n Computation successfully performed in {} seconds.".format(elapsedTime), comm)
  return 0


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':
    main()