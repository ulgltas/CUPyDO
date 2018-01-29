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
FSICouplerPath = os.path.join(os.path.dirname(__file__), '../../')
#appsPath = os.path.join(os.path.dirname(__file__), '../../apps')

sys.path.append(FSICouplerPath)
#sys.path.append(appsPath)

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
    p['nDim'] = 3
    p['tollFSI'] = 1e-6
    p['dt'] = 0.0
    p['tTot'] = 0.05
    p['nFSIIterMax'] = 4
    p['rbfRadius'] = 0.3
    p['timeIterTreshold'] = -1
    p['omegaMax'] = 1.0
    p['computationType'] = 'steady'
    p['mtfSaveAllFacs'] = False
    p['nodalLoadsType'] = 'force'
    p['nZones_SU2'] = 0
    p['withMPI'] = True
    p.update(_p)
    return p 

def main(_p, nogui):

    p = getParameters(_p)

    comm = None
    myid = 0
    numberPart = 0
    rootProcess = 0

    cupyutil.load(fileName, p['withMPI'], comm, myid, numberPart)

    if p['withMPI']:
      from mpi4py import MPI
      comm = MPI.COMM_WORLD
      myid = comm.Get_rank()
      numberPart = comm.Get_size()
    else:
      comm = None
      myid = 0
      numberPart = 1

    # --- Input parameters --- #
    CFD_file = '../../../tests/SU2_Metafor/AGARD445_Static_SU2Conf.cfg'
    CSD_file = 'AGARD445_Static_MetaforConf'

    # --- Initialize the fluid solver --- #
    import cupydoInterfaces.SU2Interface
    if comm != None:
        FluidSolver = cupydoInterfaces.SU2Interface.SU2Solver(CFD_file, p['nZones_SU2'], p['nDim'], p['computationType'], p['nodalLoadsType'], p['withMPI'], comm)
    else:
        FluidSolver = cupydoInterfaces.SU2Interface.SU2Solver(CFD_file, p['nZones_SU2'], p['nDim'], p['computationType'], p['nodalLoadsType'], p['withMPI'], 0)
    cupyutil.mpiBarrier(comm)

    # --- Initialize the solid solver --- #
    SolidSolver = None
    if myid == rootProcess:
        import cupydoInterfaces.MtfInterface
        SolidSolver = cupydoInterfaces.MtfInterface.MtfSolver(CSD_file, p['computationType'])
        SolidSolver.saveAllFacs = p['mtfSaveAllFacs']
    cupyutil.mpiBarrier(comm)

    # --- Initialize the FSI manager --- #
    manager = cupyman.Manager(FluidSolver, SolidSolver, p['nDim'], p['computationType'], comm)
    cupyutil.mpiBarrier()

    # --- Initialize the interpolator --- #
    interpolator = cupyinterp.RBFInterpolator(manager, FluidSolver, SolidSolver, p['rbfRadius'], comm)
    solverList = interpolator.getLinearSolvers()
    for ii in range(2):
        solverList[ii].setMaxNumberIterations(1000)
        solverList[ii].setPreconditioner("JACOBI")

    # --- Initialize the FSI criterion --- #
    criterion = cupycrit.DispNormCriterion(p['tollFSI'])
    cupyutil.mpiBarrier()

    # --- Initialize the FSI algorithm --- #
    algorithm = cupyalgo.AlgorithmBGSStaticRelax(manager, FluidSolver, SolidSolver, interpolator, criterion, p['nFSIIterMax'], p['dt'], p['tTot'], p['timeIterTreshold'], p['omegaMax'], comm)

    # --- Launch the FSI computation --- #
    algorithm.run()
  
    # --- Exit computation --- #
    del manager
    del criterion
    del FluidSolver
    del SolidSolver
    del interpolator
    del algorithm
    cupyutil.mpiBarrier(comm)
    return 0


# -------------------------------------------------------------------
#  Run Main Program
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
