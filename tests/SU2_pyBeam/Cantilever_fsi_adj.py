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

import os
from shutil import copyfile # Copy file from restart.pyBeam to solution.pyBeam, hardcoded filenames in pyBeam

def test_adj(res, tol):
    import numpy as np
    from cupydo.testing import CTest, CTests
    # Read results from file
    with open("AerodynamicCoeff.ascii", 'rb') as f:
        lines = f.readlines()
    resultA = np.genfromtxt(lines[-1:], delimiter=None)
    with open("SolidSolution.ascii", 'rb') as f:
        lines = f.readlines()
    resultS = np.genfromtxt(lines[-1:], delimiter=None)
    with open("sensitivity_E.dat", 'rb') as f:
        lines = f.readlines()
    resultAdj = np.genfromtxt(lines[-1:], delimiter=None)

    # Check convergence and results
    # residual for test is 1e-8
    if (res > tol):
       print("\n\n" + "FSI residual = " + str(res) + ", FSI tolerance = " + str(tol))
       raise Exception("FSI algo failed to converge!")

    tests = CTests()
    tests.add(CTest('Lift coefficient', resultA[2], -0.536, 1e-1, False)) # rel. tol. of 10%
    tests.add(CTest('Drag coefficient', resultA[3], 3.001, 1e-1, False)) # rel. tol. of 10%
    tests.add(CTest('Displacement (Tip, Y)', resultS[3], -0.000880, 1e-1, False)) # rel. tol. of 10%
    tests.add(CTest('Displacement (Tip, X)', resultS[2], 0.00381, 1e-1, False)) # rel. tol. of 10%
    tests.add(CTest('dcd/dE', resultAdj, 1.41152371e-5, 0.05, False)) # rel. tol. of 5% (was 0.3% before), based on central FD with deltaE of 1000 Pa
    tests.run()

def getAdjP():
    """Adjoint parameters"""
    import os
    filePath = os.path.abspath(os.path.dirname(__file__))
    p = {}
    # Solvers and config files
    p['fluidSolver'] = 'SU2'
    p['solidSolver'] = 'pyBeam'
    p['cfdFile'] = os.path.join(filePath, 'config_channel_adj.cfg')
    p['csdFile'] = '../../tests/SU2_pyBeam/config_cantilever.pyBeam'
    p['computation'] = 'adjoint'
    # FSI objects
    p['interpolator'] = 'RBF'
    p['interpType'] = 'conservative'
    p['algorithm'] = 'staticBGS'
    # FSI parameters
    p['compType'] = 'steady'
    p['nDim'] = 2
    p['dt'] = 0.
    p['tTot'] = 0.05
    
    p['dtSave'] = 0
    p['tol'] = 1e-8
    p['maxIt'] = 25
    p['omega'] = 1.0
    p['rbfRadius'] = .3
    p['nodalLoadsType'] = 'force'
    # Solid parameters
    p['extractors'] = [20]
    # Coupling Type

    p['mechanical'] = True
    p['thermal'] = False
    return p

def main():
    direct_case = "../../tests/SU2_pyBeam/Cantilever_fsi.py"
    exec(open(direct_case, 'r', encoding='utf8').read(), globals()) # Run direct
    copyfile('restart.pyBeam', 'solution.pyBeam')
    import cupydo.interfaces.Cupydo as cupy
    p = getAdjP() # get parameters
    cupydo = cupy.CUPyDO(p) # create fsi driver
    cupydo.run() # run fsi process
    test_adj(cupydo.algorithm.errValue, p['tol']) # check the results
    
    # eof
    print('')

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':
    main()
