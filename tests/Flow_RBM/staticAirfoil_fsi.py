#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys

filePath = os.path.abspath(os.path.dirname(sys.argv[0]))
fileName = os.path.splitext(os.path.basename(__file__))[0]
FSICouplerPath = os.path.join(os.path.dirname(__file__), '../../')

sys.path.append(FSICouplerPath)

from math import *
from optparse import OptionParser

import cupydo.utilities as cupyutil
import cupydo.manager as cupyman
import cupydo.interpolator as cupyinterp
import cupydo.criterion as cupycrit
import cupydo.algorithm as cupyalgo

from cupydo.testing import *
import numpy as np

def test(nogui, res, tol):
    # Flow constant (defined in case_fluid)
    dynP = 0.5*100 # dynamic pressure
    alpha0 = np.radians(3) # initial angle of attack
    cRef = 1

    # RBM constant (defined in case_solid)
    k = 250 # vertical stiffness
    kappa = 25 # rotational stiffness

    # Airfoil slopes (measured from flow sovler)
    cl_alpha = 6.8
    cd_alpha = 0.085
    cm_alpha = 0.27 # must be measured from flexural axis posit. hard to calibrate but of crucial importance!

    # Pitch-plunge system of equations
    A = np.array([ [k, -dynP*cRef*(cl_alpha*np.cos(alpha0) - cd_alpha*np.sin(alpha0))],
                   [0, kappa - dynP*cRef*cRef*cm_alpha] ])
    b = np.array([0, kappa*alpha0])
    x = np.linalg.solve(A, b)

    # Display the solution
    if not nogui:
        print "Ref. lift coefficient: " + str(cl_alpha*x[1])
        print "Ref. vertical displacement: " + str(x[0])
        print "Ref. new angle of attack: " + str(np.degrees(x[1]))
        print "Ref. rotational displacement : " + str(np.degrees(x[1]-alpha0))

    # Read results from file
    resultS = np.genfromtxt("NativeHistory.dat", delimiter=None, skip_header=1)
    with open("FlowHistory.dat", 'rb') as f:
        lines = f.readlines()
    resultA = np.genfromtxt(lines[-1:], delimiter=None)

    # Check convergence and results
    if (res > tol):
        print "\n\n" + "FSI residual = " + str(res) + ", FSI tolerance = " + str(tol)
        raise Exception(ccolors.ANSI_RED + "FSI algo failed to converge!" + ccolors.ANSI_RESET)
    tests = CTests()
    tests.add(CTest('Lift coefficient', resultA[2], cl_alpha*x[1], 5e-2, False)) # rel. tol. of 5%
    tests.add(CTest('Vertical displacement', -resultS[2]/cRef, x[0]/cRef, 1e-2, True)) # abs. tol. of 1% of chord
    tests.add(CTest('Rotational displacement', np.degrees(resultS[3]), np.degrees(x[1]-alpha0), 5e-1, True)) # abs. tol. of .5°
    tests.run()

def getParameters(_p):
    # --- Input parameters --- #
    p = {}
    p['nthreads'] = 1
    p['nDim'] = 2
    p['tollFSI'] = 1e-6
    p['dt'] = 0.0
    p['tTot'] = 0.0
    p['nFSIIterMax'] = 50
    p['timeIterTreshold'] = -1
    p['omegaMax'] = 1.0
    p['nbTimeToKeep'] = 0
    p['computeTangentMatrixBasedOnFirstIt'] = False
    p['computationType'] = 'steady'
    p.update(_p)
    return p

def main(_p, nogui):
    
    # --- Get FSI parameters ---#
    p = getParameters(_p)

    # --- Set up MPI and workspace --- #
    withMPI, comm, myid, numberPart = cupyutil.getMpi()
    rootProcess = 0
    cupyutil.load(fileName, withMPI, comm, myid, numberPart)
    
    # --- Input files --- #
    cfd_module = fileName[:-3] + "fluid"
    csd_file = filePath + "/" + fileName[:-3] + "solid.cfg"
    
    # --- Initialize the fluid solver --- #
    import cupydoInterfaces.FlowInterface
    fluidSolver = cupydoInterfaces.FlowInterface.Flow(cfd_module, p['nthreads'])
    
    cupyutil.mpiBarrier(comm)
    
    # --- Initialize the solid solver --- #
    solidSolver = None
    if myid == rootProcess:
        import cupydoInterfaces.RBMIntegratorInterface
        solidSolver = cupydoInterfaces.RBMIntegratorInterface.RBMIntegrator(csd_file, p['computationType'])
        
    cupyutil.mpiBarrier(comm)
        
    # --- Initialize the FSI manager --- #
    manager = cupyman.Manager(fluidSolver, solidSolver, p['nDim'], p['computationType'], comm)
    cupyutil.mpiBarrier()

    # --- Initialize the interpolator --- #
    interpolator = cupyinterp.MatchingMeshesInterpolator(manager, fluidSolver, solidSolver, comm)
    #interpolator = cupyinterp.RBFInterpolator(manager, fluidSolver, solidSolver, 1., comm)
    #interpolator = cupyinterp.TPSInterpolator(manager, fluidSolver, solidSolver, comm)
    
    # --- Initialize the FSI criterion --- #
    criterion = cupycrit.DispNormCriterion(p['tollFSI'])
    cupyutil.mpiBarrier()
    
    # --- Initialize the FSI algorithm --- #
    #algorithm = cupyalgo.AlgorithmBGSStaticRelax(manager, fluidSolver, solidSolver, interpolator, criterion, p['nFSIIterMax'], p['dt'], p['tTot'], p['timeIterTreshold'], p['omegaMax'], comm)
    #algorithm = cupyalgo.AlgorithmBGSAitkenRelax(manager, fluidSolver, solidSolver, interpolator, criterion, p['nFSIIterMax'], p['dt'], p['tTot'], p['timeIterTreshold'], p['omegaMax'], comm)
    algorithm = cupyalgo.AlgorithmIQN_ILS(manager, fluidSolver, solidSolver, interpolator, criterion, p['nFSIIterMax'], p['dt'], p['tTot'], p['timeIterTreshold'], p['omegaMax'], p['nbTimeToKeep'], p['computeTangentMatrixBasedOnFirstIt'], comm)
    
    # --- Launch the FSI computation --- #
    algorithm.run()

    # --- Check the results --- #
    test(nogui, algorithm.errValue, p['tollFSI'])

    # eof
    print ''

# -------------------------------------------------------------------
#    Run Main Program
# -------------------------------------------------------------------

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':
    
    p = {}
    
    parser=OptionParser()
    parser.add_option("--nogui", action="store_true",
                        help="Specify if we need to use the GUI", dest="nogui", default=False)
    parser.add_option("--nthreads", type="int", help="Number of threads", dest="nthreads", default=1)
    
    
    (options, args)=parser.parse_args()
    
    nogui = options.nogui
    p['nthreads'] = options.nthreads
    
    main(p, nogui)