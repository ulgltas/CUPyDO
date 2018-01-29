#! /usr/bin/env python
# -*- coding: latin-1; -*-

''' 

Copyright 2018 University of Liège

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 

'''

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