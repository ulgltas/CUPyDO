#! /usr/bin/env python3
# -*- coding: utf8 -*-

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

filePath = os.path.abspath(os.path.dirname(__file__))
fileName = os.path.splitext(os.path.basename(__file__))[0]


from math import *
from optparse import OptionParser

import cupydo.utilities as cupyutil
import cupydo.manager as cupyman
import cupydo.interpolator as cupyinterp
import cupydo.criterion as cupycrit
import cupydo.algorithm as cupyalgo

def main():

    raise Exception('Test is not working (free(): invalid pointer after Metafor sucessfull run). AFAIK, it was not tested anymore when I started with CUPyDO.\n')
    
    # --- Set up MPI --- #
    withMPI, comm, myid, numberPart = cupyutil.getMpi()
    rootProcess = 0
    
    # --- Initialize the solid solver --- #
    solid = None
    if myid == rootProcess:
        import cupydo.interfaces.Metafor as sItf
        solid = sItf.Metafor('beam', 'unsteady')
    cupyutil.mpiBarrier(comm)
    
    # --- Initialize the FSI algorithm --- #
    fsi_algo = cupyalgo.FsiSolidTestAlgorithm(solid)
    
    # --- Launch the FSI computation --- #
    fsi_algo.run()
    
    # --- Exit the solid solver --- #
    if myid == rootProcess:
        solid.exit()
    
    # --- Exit computation --- #
    cupyutil.mpiBarrier(comm)
    return 0
    
if __name__ == "__main__":
    main()
