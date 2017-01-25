import os, sys

filePath = os.path.abspath(os.path.dirname(sys.argv[0]))
fileName = os.path.splitext(os.path.basename(__file__))[0]
fsiPath = os.path.abspath('..')

sys.path.append(fsiPath)

from math import *
from optparse import OptionParser
import FSICoupler as fsi

def main():
    
    # --- Workspace set up --- #
    withMPI = False
    comm = None
    myid = 0
    numberPart = 0
    rootProcess = 0
    
    fsi.load(fileName, withMPI, comm, myid, numberPart)
    
    # --- Initialize the solid solver --- #
    solid = None
    if myid == rootProcess:
        import MtfInterface
        solid = MtfInterface.MtfSolver('beam')
    fsi.mpiBarrier(comm)
    
    # --- Initialize the FSI algorithm --- #
    fsi_algo = fsi.FsiSolidTestAlgorithm(solid)
    
    # --- Launch the FSI computation --- #
    fsi_algo.run()
    
    # --- Exit the solid solver --- #
    if myid == rootProcess:
        solid.exit()
    
    # --- Exit computation --- #
    fsi.mpiBarrier(comm)
    return 0
    
if __name__ == "__main__":
    main()