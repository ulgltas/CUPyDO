#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# original name: 

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

def test(res, tol, it):
    import numpy as np
    from cupydo.testing import CTest, CTests
    # Check convergence and results
    if (res > tol):
        print("\n\n" + "FSI residual = " + str(res) + ", FSI tolerance = " + str(tol))
        raise Exception("FSI algo failed to converge!")
    
    # Read results from file
    with open("Node_9_POS.ascii", 'rb') as f:
        lines = f.readlines()
    result_1 = np.genfromtxt(lines[-1:], delimiter=None)
    
    tests = CTests()
    tests.add(CTest('Mean nb of FSI iterations', it, 3, 1, True))
    tests.add(CTest('X-coordinate gate tip', result_1[0], 0.487164, 1e-2, False))
    tests.add(CTest('Y-coordinate gate tip', result_1[1], 0.0006501, 1e-2, False))
    tests.run()

def getFsiP():
    """Fsi parameters"""
    import os
    p = {}
    # Solvers and config files
    p['fluidSolver'] = 'Pfem'
    p['solidSolver'] = 'Metafor'
    p['cfdFile'] = 'rho1100_fluid'
    p['csdFile'] = 'rho1100_solid'
    # FSI objects
    p['interpolator'] = 'Matching'
    p['criterion'] = 'Displacements'
    p['algorithm'] = 'IQN_ILS'
    # FSI parameters
    p['compType'] = 'unsteady'
    p['computation'] = 'direct'
    p['nDim'] = 2
    p['dt'] = 0.001
    p['tTot'] = 0.05
    
    p['dtSave'] = 0
    p['tol'] = 1e-6
    p['maxIt'] = 10
    p['omega'] = 0.5
    p['firstItTgtMat'] = False
    p['nSteps'] = 6
    return p

def main():
    import cupydo.interfaces.Cupydo as cupy
    p = getFsiP() # get parameters
    cupydo = cupy.CUPyDO(p) # create fsi driver
    cupydo.run() # run fsi process
    test(cupydo.algorithm.errValue, p['tol'], cupydo.algorithm.getMeanNbOfFSIIt()) # check the results
    
    # eof
    print('')

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':
    main()
