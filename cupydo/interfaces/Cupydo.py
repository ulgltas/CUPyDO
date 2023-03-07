#! /usr/bin/env python3
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
from .. import utilities as cupyutil
from .. import manager as cupyman
from .. import interpolator as cupyinterp
from .. import criterion as cupycrit
from .. import algorithm as cupyalgo

class CUPyDO(object):
    def __init__(self, p):

        # --- Reads the solver paths --- #
        path = cupyutil.solverPath()

        # --- Set up MPI --- #
        withMPI, comm, myId, numberPart = cupyutil.getMpi()

        # --- Initialize the fluid solver --- #
        path.add(p['fluidSolver'])
        fluidSolver = self.__initFluid(p, withMPI, comm)
        path.remove(p['fluidSolver'])

        # --- Initialize the solid solver --- #
        path.add(p['solidSolver'])
        cupyutil.mpiBarrier(comm)
        solidSolver = self.__initSolid(p, myId, withMPI, comm)
        cupyutil.mpiBarrier(comm)
        path.remove(p['solidSolver'])

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
            print('Explicit simulations requested, criterion redefined to None.')
        if p['criterion'] == 'Displacements':
            criterion = cupycrit.DispNormCriterion(p['tol'])
        else:
            raise RuntimeError(p['criterion'], 'not available! (avail: "Displacements").\n')
        cupyutil.mpiBarrier()

        # --- Initialize the FSI algorithm --- #
        if p['computation'] == 'Adjoint':
            if p['algorithm'] == 'StaticBGS':
                self.algorithm = cupyalgo.AlgorithmBGSStaticRelaxAdjoint(manager, fluidSolver, solidSolver, interpolator, criterion,
                    p['maxIt'], p['dt'], p['tTot'], p['dtSave'], p['omega'], comm)
            else:
                raise RuntimeError(p['algorithm'], 'not available in adjoint calculations! (avail: "StaticBGS").\n')

        else:
            if p['algorithm'] == 'Explicit':
                self.algorithm = cupyalgo.AlgorithmExplicit(manager, fluidSolver, solidSolver, interpolator,
                    p['dt'], p['tTot'], p['dtSave'], comm)
            elif p['algorithm'] == 'StaticBGS':
                self.algorithm = cupyalgo.AlgorithmBGSStaticRelax(manager, fluidSolver, solidSolver, interpolator, criterion,
                    p['maxIt'], p['dt'], p['tTot'], p['dtSave'], p['omega'], comm)
            elif p['algorithm'] == 'AitkenBGS':
                self.algorithm = cupyalgo.AlgorithmBGSAitkenRelax(manager, fluidSolver, solidSolver, interpolator, criterion,
                    p['maxIt'], p['dt'], p['tTot'], p['dtSave'], p['omega'], comm)
            elif p ['algorithm'] == 'IQN_ILS':
                self.algorithm = cupyalgo.AlgorithmIQN_ILS(manager, fluidSolver, solidSolver, interpolator, criterion,
                    p['maxIt'], p['dt'], p['tTot'], p['dtSave'], p['omega'], p['nSteps'], p['firstItTgtMat'], comm)
            elif p ['algorithm'] == 'IQN_MVJ':
                self.algorithm = cupyalgo.AlgorithmIQN_MVJ(manager, fluidSolver, solidSolver, interpolator, criterion,
                    p['maxIt'], p['dt'], p['tTot'], p['dtSave'], p['omega'], p['firstItTgtMat'], comm)
            else:
                raise RuntimeError(p['algorithm'], 'not available! (avail: "Explicit", "StaticBGS", "AitkenBGS", "IQN_ILS" or "IQN_MVJ").\n')
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
        if p['computation'] == 'Adjoint':
            if p['fluidSolver'] == 'SU2':
                from . import SU2 as fItf
                if comm != None:
                    fluidSolver = fItf.SU2Adjoint(p['cfdFile'], p['nDim'], p['compType'], p['nodalLoadsType'], withMPI, comm)
                else:
                    fluidSolver = fItf.SU2Adjoint(p['cfdFile'], p['nDim'], p['compType'], p['nodalLoadsType'], withMPI, 0)
            else:
                raise RuntimeError('Adjoint interface for', p['fluidSolver'], 'not found!\n')
        else:
            if p['fluidSolver'] == 'SU2':
                from . import SU2 as fItf
                if comm != None:
                    fluidSolver = fItf.SU2(p['cfdFile'], p['nDim'], p['compType'], p['nodalLoadsType'], withMPI, comm)
                else:
                    fluidSolver = fItf.SU2(p['cfdFile'], p['nDim'], p['compType'], p['nodalLoadsType'], withMPI, 0)
            elif p['fluidSolver'] == 'Pfem':
                from . import Pfem as fItf
                fluidSolver = fItf.Pfem(p['cfdFile'], args.n, args.nogui, p['dt'])
            elif p['fluidSolver'] == 'DART':
                from dart.api.cupydo import Dart
                fluidSolver = Dart(p['cfdFile'], args.n)
            elif p['fluidSolver'] == 'VLM':
                from . import VLM as fItf
                fluidSolver = fItf.VLMSolver(p['cfdFile'])
            elif p['fluidSolver'] == 'Pfem3D':
                from . import Pfem3D as fItf
                fluidSolver = fItf.Pfem3D(p)
            else:
                raise RuntimeError('Interface for', p['fluidSolver'], 'not found!\n')
        return fluidSolver

    def __initSolid(self, p, myId, withMPI, comm):
        """Initialize solid solver interface
        Adrien Crovato
        """
        solidSolver = None
        # IMPORTANT! only master can instantiate the solid solver except SU2Solid.
        if p['computation'] == 'Adjoint': # Adjoint calculations only support SU2Solid
            if p['solidSolver'] == 'SU2':
                from . import SU2Solid as sItf
                if comm != None:
                    solidSolver = sItf.SU2SolidAdjoint(p['csdFile'], p['nDim'], p['compType'], p['nodalLoadsType'], p['extractors'], p['surfaceFilename'], p['surfaceExtension'], withMPI, comm)
                else:
                    solidSolver = sItf.SU2SolidAdjoint(p['csdFile'], p['nDim'], p['compType'], p['nodalLoadsType'], p['extractors'], p['surfaceFilename'], p['surfaceExtension'], withMPI, 0)
            elif myId == 0 and p['solidSolver'] == 'pyBeam':
                from . import Beam as sItf
                solidSolver = sItf.pyBeamAdjointSolver(p['csdFile'], p['nDim'], p['compType'], p['nodalLoadsType'], p['extractors'])
            elif myId == 0:
                raise RuntimeError('Adjoint interface for', p['solidSolver'], 'not found!\n')
        else:
            if myId == 0 and p['solidSolver'] == 'Metafor':
                from . import Metafor as sItf
                solidSolver = sItf.Metafor(p)
            elif myId == 0 and p['solidSolver'] == 'RBMI':
                from . import RBMI as sItf
                solidSolver = sItf.RBMI(p['csdFile'], p['compType'])
            elif myId == 0 and p['solidSolver'] == 'Modal':
                from . import Modal as sItf
                solidSolver = sItf.Modal(p['csdFile'], p['compType'])
            elif myId == 0 and p['solidSolver'] == 'pyBeam':
                from . import Beam as sItf
                solidSolver = sItf.pyBeamSolver(p['csdFile'], p['nDim'], p['compType'], p['nodalLoadsType'], p['extractors'])
            elif p['solidSolver'] == 'SU2':
                from . import SU2Solid as sItf
                if comm != None:
                    solidSolver = sItf.SU2SolidSolver(p['csdFile'], p['nDim'], p['compType'], p['nodalLoadsType'], p['extractors'], p['surfaceFilename'], p['surfaceExtension'], withMPI, comm)
                else:
                    solidSolver = sItf.SU2SolidSolver(p['csdFile'], p['nDim'], p['compType'], p['nodalLoadsType'], p['extractors'], p['surfaceFilename'], p['surfaceExtension'], withMPI, 0)
            elif p['solidSolver'] == 'GetDP':
                from . import GetDP as sItf
                raise RuntimeError('GetDP interface not up-to-date!')
            elif myId == 0:
                raise RuntimeError('Interface for', p['solidSolver'], 'not found!\n')
        return solidSolver

# ---------------------------------------------------------------------- #
# Sample parameters list

# Solvers
# - p['fluidSolver'], fluid solvers available: SU2, Pfem, DART, VLM, Pfem3D
# - p['solidSolver'], solid solvers available: Metafor, RBMI, Modal, GetDP, SU2
# Configuration files
# - p['cfdFile'], path to fluid cfg file 
# - p['csdFile'], path to solid cfg file

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
# - p['dtSave'], time between each result save()
# needed by BGS
# - p['tol'], tolerance on displacements
# - p['maxIt'], maximum number of iterations
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
# - p['extractors'], SU2Solid, pyBeam
# - p['surfaceFilename'], SU2Solid, filename of the output surface file, as defined in the CFG file
# - p['surfaceExtension'], SU2Solid, extension of the output surface file: vtu, vtk, dat
# ---------------------------------------------------------------------- #
