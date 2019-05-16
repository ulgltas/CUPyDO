#! /usr/bin/env python
# -*- coding: utf-8 -*-

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

cupydo.py
Drive (create and run) an FSI computation
Authors: Adrien CROVATO

'''

import cupydo.utilities as cupyutil
import cupydo.manager as cupyman
import cupydo.interpolator as cupyinterp
import cupydo.criterion as cupycrit
import cupydo.algorithm as cupyalgo

class CUPyDO():
    def __init__(self, p):
        # --- Set up MPI --- #
        withMPI, comm, myId, numberPart = cupyutil.getMpi()
        rootProcess = 0

        # --- Initialize the fluid and solid solvers --- #
        fluidSolver = self.__initFluid(p, withMPI, comm)
        cupyutil.mpiBarrier(comm)
        solidSolver = self.__initSolid(p, myId)
        cupyutil.mpiBarrier(comm)

        # --- Initialize the FSI manager --- #
        manager = cupyman.Manager(fluidSolver, solidSolver, p['nDim'], p['compType'], comm)
        cupyutil.mpiBarrier()

        # --- Initialize the interpolator --- #
        if p['interpolator'] == 'Matching':
            interpolator = cupyinterp.MatchingMeshesInterpolator(manager, fluidSolver, solidSolver, comm)
        elif p['interpolator'] == 'RBF':
            interpolator = cupyinterp.RBFInterpolator(manager, fluidSolver, solidSolver, p['rbfRadius'], comm)
        elif p['interpolator'] == 'TPS':
            interpolator = cupyinterp.TPSInterpolator(manager, fluidSolver, solidSolver, comm)
        else:
            raise RuntimeError(p['interpolator'], 'not available! (avail: "Matching", "RBF" or "TPS").\n')
        # if petsc is used, then some options can be set
        if withMPI and 'interpOpts' in p:
            for linSolver in interpolator.getLinearSolvers():
                linSolver.setMaxNumberIterations(p['interpOpts'][0])
                linSolver.setPreconditioner(p['interpOpts'][1])

        # --- Initialize the FSI criterion --- #
        if p['algorithm'] == 'Explicit':
            print 'Explicit simulations requested, criterion redefined to None.'
        if p['criterion'] == 'Displacements':
            criterion = cupycrit.DispNormCriterion(p['tol'])
        else:
            raise RuntimeError(p['criterion'], 'not available! (avail: "Displacements").\n')
        cupyutil.mpiBarrier()

        # --- Initialize the FSI algorithm --- #
        if p['algorithm'] == 'Explicit':
            self.algorithm = cupyalgo.AlgorithmExplicit(manager, fluidSolver, solidSolver, interpolator,
                p['dt'], p['tTot'], p['timeItTresh'], comm)
        elif p['algorithm'] == 'StaticBGS':
            self.algorithm = cupyalgo.AlgorithmBGSStaticRelax(manager, fluidSolver, solidSolver, interpolator, criterion,
                p['maxIt'], p['dt'], p['tTot'], p['timeItTresh'], p['omega'], comm)
        elif p['algorithm'] == 'AitkenBGS':
            self.algorithm = cupyalgo.AlgorithmBGSAitkenRelax(manager, fluidSolver, solidSolver, interpolator, criterion,
                p['maxIt'], p['dt'], p['tTot'], p['timeItTresh'], p['omega'], comm)
        elif p ['algorithm'] == 'IQN_ILS':
            self.algorithm = cupyalgo.AlgorithmIQN_ILS(manager, fluidSolver, solidSolver, interpolator, criterion,
                p['maxIt'], p['dt'], p['tTot'], p['timeItTresh'], p['omega'], p['nSteps'], p['firstItTgtMat'], comm)
        cupyutil.mpiBarrier()

    def run(self):
        """Launch the FSI computation
        Adrien Crovato
        """
        self.algorithm.run()

    def __initFluid(self, p, withMPI, comm):
        """Initialize fluid solver interface
        Adrien Crovato
        """
        args = cupyutil.parseArgs()
        if p['fluidSolver'] == 'SU2':
            import cupydo.interfaces.SU2 as fItf
            if comm != None:
                fluidSolver = fItf.SU2(p['cfdFile'], p['nDim'], p['compType'], p['nodalLoadsType'], withMPI, comm)
            else:
                fluidSolver = fItf.SU2(p['cfdFile'], p['nDim'], p['compType'], p['nodalLoadsType'], withMPI, 0)
        elif p['fluidSolver'] == 'Pfem':
            import cupydo.interfaces.Pfem as fItf
            fluidSolver = fItf.Pfem(p['cfdFile'], args.n, args.nogui, p['dt'])
        elif p['fluidSolver'] == 'Flow':
            import cupydo.interfaces.Flow as fItf
            fluidSolver = fItf.Flow(p['cfdFile'], args.n)
        else:
            raise RuntimeError('Interface for', p['fluidSolver'], 'not found!\n')
        return fluidSolver

    def __initSolid(self, p, myId):
        """Initialize fluid solver interface
        Adrien Crovato
        """
        solidSolver = None
        if myId == 0: # only master can instantiate the solid solver
            if p['solidSolver'] == 'Metafor':
                import cupydo.interfaces.Metafor as sItf
                solidSolver = sItf.Metafor(p['csdFile'], p['compType'])
            elif p['solidSolver'] == 'RBMI':
                import cupydo.interfaces.RBMI as sItf
                solidSolver = sItf.RBMI(p['csdFile'], p['compType'])
            elif p['solidSolver'] == 'Modal':
                import cupydo.interfaces.Modal as sItf
                solidSolver = sItf.Modal(p['csdFile'], p['compType'])
            elif p['solidSolver'] == 'GetDP':
                import cupydo.interfaces.GetDP as sItf
                raise RuntimeError('GetDP interface not up-to-date!\n')
            else:
                raise RuntimeError('Interface for', p['solidSolver'], 'not found!\n')
        return solidSolver

# ---------------------------------------------------------------------- #
# Sample parameters list

# Solvers
# - p['fluidSolver'], fluid solvers available: SU2, Pfem, Flow
# - p['solidSolver'], solid solvers available: Metafor, RBMI, Modal, GetDP
# Configuration files
# - p['cfdFile'], path to fluid cfg file 
# - p['csdFile'], path to solid cfg file'

# FSI objects
# - p['interpolator'], interpolator type available: Matching, RBF, TPS
# - p['criterion'], convergence criterion available: Displacements
# - p['algorithm'], FSI algorithms available: Explicit, StaticBGS, AitkenBGS, IQN_ILS

# FSI parameters
# needed by all algos
# - p['compType'], steady or unsteady
# - p['nDim'], dimension, 2 or 3
# - p['dt'], time steps
# - p['tTot'], total time
# - p['timeItTresh'], ??????
# needed by BGS
# - p['tol'], tolerance on displacements
# - p['maxIt'], maximu number of iterations
# - p['omega'], relaxation parameter
# needed by IQN-ILS
# - p['firstItTgtMat'], compute the Tangent matrix based on first iteration (True or False)
# - p['nSteps'], number of time steps to keep
# needed by RBF interpolator
# - p['rbfRadius'], radius of interpolation for RBF
# optional for RBF/TPS interpolators
# - p['interpOpts'], optional options for interpolator, [0] = max number of iterations, [1] = preconditionner type

# Solver parameters that should be moved to solver cfg files and handled by the solver interface
# - p['nodalLoadsType'], SU2
# ---------------------------------------------------------------------- #
