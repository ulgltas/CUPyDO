#!/usr/bin/env python

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

from optparse import OptionParser	# use a parser for configuration
from math import *
import numpy as np

import cupydo.utilities as cupyutil
import cupydo.manager as cupyman
import cupydo.interpolator as cupyinterp

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():

    cfd_file = 'SU2Conf.cfg'
    csd_file = 'surface_wing_modes.csv'
    initialModalDisp = np.array([5.6023, 0., 0., 0., 0., 0.])
    initialModalVel = np.zeros(6, dtype=float)
    modalMass = np.diag([1., 1., 1., 1., 1., 1.])
    modalDamping = np.zeros((6, 6), dtype=float)
    modalStiffness = np.diag([1.817819061E+02, 1.178750000E+03, 1.709068970E+03, 5.454931152E+03, 8.660587891E+03, 9.815511719E+03])
    nzone = 1
    nDim = 3
    computationType = 'unsteady'
    nodalLoadsType = 'force'
    withMPI = True
    rootProcess = 0
    deltaT = 6.226663e-03
    totTime = 1.497512e+00

    if withMPI == True:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        myid = comm.Get_rank()
        numberPart = comm.Get_size()
    else:
        comm = 0 
        myid = 0
        numberPart = 1

    # --- Initialize the fluid solver --- #
    import cupydoInterfaces.SU2Interface
    if comm != None:
        fluidSolver = cupydoInterfaces.SU2Interface.SU2Solver(cfd_file, nzone, nDim, computationType, nodalLoadsType, withMPI, comm)
    else:
        fluidSolver = cupydoInterfaces.SU2Interface.SU2Solver(cfd_file, nzone, nDim, computationType, nodalLoadsType, withMPI, 0)

    # --- Initialize modal interpreter --- #
    solidSolver = None
    if myid == rootProcess:
        import ModalSolverInterface
        solidSolver = ModalSolverInterface.modalInterpreter(csd_file, initialModalDisp, initialModalVel, modalMass, modalDamping, modalStiffness)

    cupyutil.mpiBarrier(comm)

    # --- Initialize the manager --- #
    manager = cupyman.Manager(fluidSolver, solidSolver, nDim, computationType, comm)
    cupyutil.mpiBarrier()

    # --- Initialize the interpolator --- #
    interpolator = cupyinterp.MatchingMeshesInterpolator(manager, fluidSolver, solidSolver, comm)

    # --- Run the computation --- #
    timeIter = 0
    nTimeIter = int((totTime/deltaT)-1)
    time = timeIter*deltaT

    while(timeIter <= nTimeIter):
        fluidSolver.preprocessTimeIter(timeIter)
        if myid in manager.getSolidSolverProcessors():
            solidSolver.run(time-deltaT, time)
        interpolator.getDisplacementFromSolidSolver()
        interpolator.interpolateSolidDisplacementOnFluidMesh()
        interpolator.setDisplacementToFluidSolver(time)
        fluidSolver.meshUpdate(timeIter)
        fluidSolver.run(time-deltaT, time)
        
        interpolator.getLoadsFromFluidSolver()
        interpolator.interpolateFluidLoadsOnSolidMesh()
        interpolator.setLoadsToSolidSolver(time)        

        fluidSolver.update(deltaT)
        
        stopComp = fluidSolver.save(timeIter)
        if (stopComp == True):
            break
        timeIter += 1
        time += deltaT

    # --- Exit fluid solver --- #
    fluidSolver.exit()

    # --- Exit solid solver --- #
    if myid == rootProcess:
        solidSolver.exit()


    del fluidSolver
    del solidSolver
    del interpolator
    del manager

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()  
