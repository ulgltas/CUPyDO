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
        manager = cupyman.Manager(fluidSolver, solidSolver, p, comm)
        cupyutil.mpiBarrier()

        # --- Initialize the interpolator --- #
        if p['interpolator'] == 'matching':
            interpolator = cupyinterp.MatchingMeshesInterpolator(manager, fluidSolver, solidSolver, p, comm)
        elif p['interpolator'] == 'RBF':
            if p['interpType'] == 'consistent':
                interpolator = cupyinterp.ConsistentRBFInterpolator(manager, fluidSolver, solidSolver, p, comm)
            elif p['interpType'] == 'conservative':
                interpolator = cupyinterp.ConservativeRBFInterpolator(manager, fluidSolver, solidSolver, p, comm)
            else:
                raise RuntimeError(p['interpType'], 'not available! (avail: conservative or consistent).\n')
        elif p['interpolator'] == 'TPS':
            if p['interpType'] == 'consistent':
                interpolator = cupyinterp.ConsistentTPSInterpolator(manager, fluidSolver, solidSolver, p, comm)
            elif p['interpType'] == 'conservative':
                interpolator = cupyinterp.ConservativeTPSInterpolator(manager, fluidSolver, solidSolver, p, comm)
            else:
                raise RuntimeError(p['interpType'], 'not available! (avail: conservative or consistent).\n')
        else:
            raise RuntimeError(p['interpolator'], 'not available! (avail: matching, RBF or TPS).\n')
        # if petsc is used, then some options can be set
        if withMPI and 'interpOpts' in p:
            for linSolver in interpolator.getLinearSolvers():
                linSolver.setMaxNumberIterations(p['interpOpts'][0])
                linSolver.setPreconditioner(p['interpOpts'][1])

        # --- Initialize the FSI criterion --- #

        criterion = cupycrit.NormCriterion(p)
        cupyutil.mpiBarrier()

        # --- Initialize the FSI algorithm --- #
        if p['computation'] == 'adjoint':
            if p['algorithm'] == 'staticBGS':
                self.algorithm = cupyalgo.AlgorithmBGSStaticRelaxAdjoint(manager, fluidSolver, solidSolver, interpolator, criterion, p, comm)
            else:
                raise RuntimeError(p['algorithm'], 'not available in adjoint calculations! (avail: staticBGS).\n')

        else:
            if p['algorithm'] == 'explicit':
                self.algorithm = cupyalgo.AlgorithmExplicit(manager, fluidSolver, solidSolver, interpolator, p, comm)
            elif p['algorithm'] == 'staticBGS':
                self.algorithm = cupyalgo.AlgorithmBGSStaticRelax(manager, fluidSolver, solidSolver, interpolator, criterion, p, comm)
            elif p['algorithm'] == 'aitkenBGS':
                self.algorithm = cupyalgo.AlgorithmBGSAitkenRelax(manager, fluidSolver, solidSolver, interpolator, criterion, p, comm)
            elif p ['algorithm'] == 'IQN_ILS':
                self.algorithm = cupyalgo.AlgorithmIQN_ILS(manager, fluidSolver, solidSolver, interpolator, criterion, p, comm)
            elif p ['algorithm'] == 'IQN_MVJ':
                self.algorithm = cupyalgo.AlgorithmIQN_MVJ(manager, fluidSolver, solidSolver, interpolator, criterion, p, comm)
            else:
                raise RuntimeError(p['algorithm'], 'not available! (avail: explicit, staticBGS, aitkenBGS, IQN_ILS or IQN_MVJ).\n')
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
        if p['computation'] == 'adjoint':
            if p['fluidSolver'] == 'SU2':
                from . import SU2 as fItf
                if comm != None:
                    fluidSolver = fItf.SU2Adjoint(p, withMPI, comm)
                else:
                    fluidSolver = fItf.SU2Adjoint(p, withMPI, 0)
            else:
                raise RuntimeError('Adjoint interface for', p['fluidSolver'], 'not found!\n')
        else:
            if p['fluidSolver'] == 'SU2':
                from . import SU2 as fItf
                if comm != None:
                    fluidSolver = fItf.SU2(p, withMPI, comm)
                else:
                    fluidSolver = fItf.SU2(p, withMPI, 0)
            elif p['fluidSolver'] == 'Pfem':
                from . import Pfem as fItf
                fluidSolver = fItf.Pfem(p, args.k)
            elif p['fluidSolver'] == 'DART':
                from dart.api.cupydo import Dart
                fluidSolver = Dart(p['cfdFile'], args.k)
            elif p['fluidSolver'] == 'VLM':
                from . import VLM as fItf
                fluidSolver = fItf.VLMSolver(p)
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
        if p['computation'] == 'adjoint': # Adjoint calculations only support SU2Solid
            if p['solidSolver'] == 'SU2':
                from . import SU2Solid as sItf
                if comm != None:
                    solidSolver = sItf.SU2SolidAdjoint(p, withMPI, comm)
                else:
                    solidSolver = sItf.SU2SolidAdjoint(p, withMPI, 0)
            elif myId == 0 and p['solidSolver'] == 'pyBeam':
                from . import Beam as sItf
                solidSolver = sItf.pyBeamAdjointSolver(p)
            elif myId == 0:
                raise RuntimeError('Adjoint interface for', p['solidSolver'], 'not found!\n')
        else:
            if myId == 0 and p['solidSolver'] == 'Metafor':
                from . import Metafor as sItf
                solidSolver = sItf.Metafor(p)
            elif myId == 0 and p['solidSolver'] == 'RBMI':
                from . import RBMI as sItf
                solidSolver = sItf.RBMI(p)
            elif myId == 0 and p['solidSolver'] == 'Modal':
                from . import Modal as sItf
                solidSolver = sItf.Modal(p)
            elif myId == 0 and p['solidSolver'] == 'pyBeam':
                from . import Beam as sItf
                solidSolver = sItf.pyBeamSolver(p)
            elif p['solidSolver'] == 'SU2':
                from . import SU2Solid as sItf
                if comm != None:
                    solidSolver = sItf.SU2SolidSolver(p, withMPI, comm)
                else:
                    solidSolver = sItf.SU2SolidSolver(p, withMPI, 0)
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
# - p['interpType'], loading type available: Conservative, Consistent
# - p['criterion'], convergence criterion available: Displacements
# - p['algorithm'], FSI algorithms available: Explicit, StaticBGS, AitkenBGS, IQN_ILS

# FSI parameters
# needed by all algos
# - p['regime'], steady or unsteady
# - p['computation'], direct or adjoint
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
# - p['extractors'], SU2Solid, pyBeam
# - p['surfaceFilename'], SU2Solid, filename of the output surface file, as defined in the CFG file
# - p['surfaceExtension'], SU2Solid, extension of the output surface file: vtu, vtk, dat
# ---------------------------------------------------------------------- #
