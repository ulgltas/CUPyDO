#! /usr/bin/env python
# -*- coding: utf8 -*-

'''

Copyright 2025 Université de Liège

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.'''

import cupydo.interfaces.Cupydo as cupy
import numpy as np
from cupydo.utilities import mpiPrint, getMpi, mpiAllReduce, mpiBarrier
import os

omegaHB = 63.0
haveMPI, comm, myid, numberPart = getMpi()

def test(res, tol):
    global obj_fun
    import numpy as np
    from cupydo.testing import CTest, CTests
    # Read results from file
    with open("NativeHistoryFSI.dat", 'rb') as f:
        lines = f.readlines()
    resultRBM = np.genfromtxt(lines[-3:], delimiter=None) # Last three time instances of the RBMI HB solution

    # Check convergence and results
    if (res > tol):
        print("\n\n" + "FSI residual = " + str(res) + ", FSI tolerance = " + str(tol))
        raise Exception("FSI algo failed to converge!")
    tests = CTests()
    tests.add(CTest('Mean drag', obj_fun, 0.0337, 1e-1, False)) # rel. tol. of 10%
    tests.add(CTest('First instance pitch displacement', resultRBM[0, 3], 0.0792, 1e-1, False)) # rel. tol. of 10%
    tests.run()

def getFsiP():
    """Fsi parameters"""
    import os
    filePath = os.path.abspath(os.path.dirname(__file__))
    p = {}
    # Solvers and config files

    p['fluidSolver'] = 'SU2'
    p['solidSolver'] = 'RBMI'
    p['cfdFile'] = os.path.join(filePath, 'SU2Conf_V_080_HB_1.cfg')
    p['csdFile'] = os.path.join(filePath, 'PitchPlungeAirfoil_HB_1_RBMConf.cfg')
    p['computation'] = 'direct'

    # FSI objects

    p['interpolator'] = 'matching'
    p['interpType'] = 'conservative'
    p['algorithm'] = 'staticBGS'

    # FSI parameters

    p['regime'] = 'harmonic'
    p['omegaHB'] = 63.0
    p['criterion'] = 'relative'
    p['nInst'] = 3
    p['nDim'] = 2
    p['dt'] = 0.001
    p['tTot'] = 0.005
    p['dtSave'] = 0
    p['maxIt'] = 20
    p['omega'] = 0.7

    # Coupling Type

    p['mechanical'] = True
    p['mechanicalTol'] = 1e-3
    p['thermal'] = False
    p['omegaTol'] = 1e-3
    return p

def getFsiAdjP():
    """Fsi adjoint parameters"""
    import os
    filePath = os.path.abspath(os.path.dirname(__file__))
    p = {}
    # Solvers and config files

    p['fluidSolver'] = 'SU2'
    p['solidSolver'] = 'RBMI'
    p['cfdFile'] = os.path.join(filePath, 'SU2Conf_Adjoint_Drag_HB_1.cfg')
    p['csdFile'] = os.path.join(filePath, 'AdjointDrag_HB_1_RBMConf.cfg')
    p['computation'] = 'adjoint'

    # FSI objects

    p['interpolator'] = 'matching'
    p['interpType'] = 'conservative'
    p['algorithm'] = 'staticBGS'

    # FSI parameters

    p['regime'] = 'harmonic'
    p['omegaHB'] = 63.0
    p['criterion'] = 'relative'
    p['nInst'] = 3
    p['nDim'] = 2
    p['dt'] = 0.001
    p['tTot'] = 0.005
    p['dtSave'] = 0
    p['maxIt'] = 20
    p['omega'] = 0.7

    # Coupling Type

    p['mechanical'] = True
    p['mechanicalTol'] = 5e-4
    p['thermal'] = False
    p['omegaTol'] = np.inf
    return p

def runDirect(x):
    global omegaHB
    
    mpiPrint("Running direct simulation with design variable vector: {}".format(x), comm)
    if myid == 0:
        f = open("optimisation.txt", "a")
        print("Running direct simulation with design variable vector: {}".format(np.array2string(x, max_line_width=255)), file=f)
        f.close()
    
    pDir = getFsiP() # get parameters
    pDir['omegaHB'] = omegaHB # Update frequency
    cupydoDir = cupy.CUPyDO(pDir) # create fsi driver
    if myid == 0:
        cupydoDir.algorithm.SolidSolver.applyDesignVariables(x)
    cupydoDir.run()
    if myid == 0:
        omega = cupydoDir.algorithm.omegaHB
        value = cupydoDir.algorithm.objectiveFunction
        omegaHB = mpiAllReduce(mpiComm=comm, value=omega)
        mpiBarrier(comm)
        valueOF = mpiAllReduce(mpiComm=comm, value=value)
        mpiBarrier(comm)
    else:
        omegaHB = mpiAllReduce(mpiComm=comm, value=0.)
        mpiBarrier(comm)
        valueOF = mpiAllReduce(mpiComm=comm, value=0.)
        mpiBarrier(comm)
    mpiPrint("Objective function value: {}, frequency value: {}".format(valueOF, omegaHB), comm)
    if myid == 0:
        f = open("optimisation.txt", "a")
        print("Objective function value: {}, frequency value: {}".format(valueOF, omegaHB), file=f)
        f.close()
    return valueOF

def runAdjoint(x):
    global omegaHB, resAdj
    mpiPrint("Running adjoint simulation with design variable vector: {}".format(x), comm)
    if myid == 0:
        f = open("optimisation.txt", "a")
        print("Running adjoint simulation with design variable vector: {}".format(np.array2string(x, max_line_width=255)), file=f)
        f.close()
    pAdj = getFsiAdjP() # get parameters
    pAdj['omegaHB'] = omegaHB # Update frequency
    cupydoAdj = cupy.CUPyDO(pAdj) # create fsi driver
    if myid == 0:
        cupydoAdj.algorithm.SolidSolver.applyDesignVariables(x)
    cupydoAdj.run()
    resAdj = cupydoAdj.algorithm.criterion.epsilon

    if myid == 0:
        sendGradients = cupydoAdj.algorithm.gradients
    else:
        sendGradients = np.zeros((24,))

    if comm != None:
        from mpi4py import MPI
        valueGradients = np.zeros((24,))
        comm.Allreduce(sendGradients, valueGradients, MPI.SUM)
        mpiBarrier(comm)
    else:
        valueGradients = sendGradients
    mpiPrint("Gradients value: {}, frequency value: {}".format(np.array2string(valueGradients, separator=','), omegaHB), comm)
    if myid == 0:
        f = open("optimisation.txt", "a")
        print("Gradients value: {}, frequency value: {}".format(np.array2string(valueGradients, separator=',', max_line_width=255), omegaHB), file=f)
        f.close()
    return valueGradients


def main():
    global resAdj, obj_fun
    resAdj = 0.0
    obj_fun = 0.0
    alpha0 = np.zeros((24,))
    obj_fun = runDirect(alpha0)
    gradients = runAdjoint(alpha0)

    p = getFsiAdjP()
    test(resAdj, p['mechanicalTol'])
    # eof
    print('')

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':
    main()
