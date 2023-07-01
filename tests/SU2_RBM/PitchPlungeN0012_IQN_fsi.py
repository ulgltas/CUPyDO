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

    # Check convergence and results
    if (res > tol):
        print("\n\n" + "FSI residual = " + str(res) + ", FSI tolerance = " + str(tol))
        raise Exception("FSI algo failed to converge!")
    tests = CTests()
    tests.add(CTest('Lift coefficient', resultA[2], 0.345479, 1e-1, False))
    tests.add(CTest('Drag coefficient', resultA[3], 0.022232, 1e-1, False))
    tests.run()

def getFsiP():
    """Fsi parameters"""
    import os
    filePath = os.path.abspath(os.path.dirname(__file__))
    p = {}
    # Solvers and config files
    p['fluidSolver'] = 'SU2'
    p['solidSolver'] = 'RBMI'
    p['cfdFile'] = os.path.join(filePath, 'PitchPlungeN0012_IQN_fsi_SU2Conf.cfg')
    p['csdFile'] = os.path.join(filePath, 'PitchPlungeN0012_IQN_fsi_RBMConf.cfg')
    # FSI objects
    p['interpolator'] = 'RBF'
    p['criterion'] = 'displacement'
    p['algorithm'] = 'IQN_ILS'
    # FSI parameters
    p['compType'] = 'unsteady'
    p['computation'] = 'direct'
    p['nDim'] = 2
    p['dt'] = 0.0005
    p['tTot'] = 0.005
    
    p['dtSave'] = 0
    p['tol'] = 1e-6
    p['maxIt'] = 25
    p['omega'] = 1.0
    p['nSteps'] = 0
    p['firstItTgtMat'] = False
    p['rbfRadius'] = .1
    p['nodalLoadsType'] = 'force'
    return p

def main():
    import cupydo.interfaces.Cupydo as cupy
    p = getFsiP() # get parameters
    cupydo = cupy.CUPyDO(p) # create fsi driver
    cupydo.run() # run fsi process
    test(cupydo.algorithm.errValue, p['tol']) # check the results
    
    # eof
    print('')

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':
    main()