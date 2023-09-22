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


def test(res, tol):
    import numpy as np
    from cupydo.testing import CTest, CTests
    # Read results from file
    with open("AerodynamicCoeff.ascii", 'rb') as f:
        lines = f.readlines()
    resultA = np.genfromtxt(lines[-1:], delimiter=None)
    with open("db_Field(TZ,RE)_GROUP_ID_180.ascii", 'rb') as f:
        lines = f.readlines()
    resultS1 = np.genfromtxt(lines[-1:], delimiter=None)
    with open("db_Field(TZ,RE)_GROUP_ID_181.ascii", 'rb') as f:
        lines = f.readlines()
    resultS2 = np.genfromtxt(lines[-1:], delimiter=None)

    # Check convergence and results
    # residual for test is 0.0017
    #if (res > tol):
    #    print "\n\n" + "FSI residual = " + str(res) + ", FSI tolerance = " + str(tol)
    #    raise Exception("FSI algo failed to converge!")
    tests = CTests()
    tests.add(CTest('Lift coefficient', resultA[2], 0.054409, 1e-1, False))
    tests.add(CTest('Drag coefficient', resultA[3], 0.000042, 1e-1, False))
    tests.add(CTest('Displacement (180, TZ)', resultS1[-1], 0.012873, 1e-1, False))
    tests.add(CTest('Displacement (181, TZ)', resultS2[-1], 0.014388, 1e-1, False))
    tests.run()

def getFsiP():
    """Fsi parameters"""
    import os
    filePath = os.path.abspath(os.path.dirname(__file__))
    p = {}
    # Solvers and config files
    p['fluidSolver'] = 'SU2'
    p['solidSolver'] = 'Metafor'
    p['cfdFile'] = os.path.join(filePath, 'AGARD445_Static_SU2Conf.cfg')
    p['csdFile'] = 'AGARD445_Static_MetaforConf'
    # FSI objects
    p['interpolator'] = 'TPS'
    p['interpType'] = 'conservative'
    p['algorithm'] = 'staticBGS'
    # FSI parameters
    p['regime'] = 'steady'
    p['computation'] = 'direct'
    p['nDim'] = 3
    p['dt'] = 0.
    p['tTot'] = 0.05
    
    p['dtSave'] = 0
    p['tol'] = 1e-6
    p['maxIt'] = 4
    p['omega'] = 1.0
    p['interpOpts'] = [1000, 'JACOBI']
    p['nodalLoadsType'] = 'force'
    # Coupling Type

    p['mechanical'] = True
    p['thermal'] = False
    return p

def main():
    import cupydo.interfaces.Cupydo as cupy
    p = getFsiP() # get parameters
    cupydo = cupy.CUPyDO(p) # create fsi driver
    cupydo.run() # run fsi process
    test(cupydo.algorithm.criterion.epsilon, p['tol']) # check the results
    
    # eof
    print('')

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':
    main()