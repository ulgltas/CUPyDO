#! /usr/bin/env python
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
from __future__ import print_function

def test(res, tol):
    import numpy as np
    from cupydo.testing import CTest, CTests, ccolors
    # Read results from file
    with open("AerodynamicCoeff.ascii", 'rb') as f:
        lines = f.readlines()
    resultA = np.genfromtxt(lines[-1:], delimiter=None)
    with open("SolidSolution.ascii", 'rb') as f:
        lines = f.readlines()
    resultS = np.genfromtxt(lines[-1:], delimiter=None)

    # Check convergence and results
    # residual for test is 0.002
    #if (res > tol):
    #    print "\n\n" + "FSI residual = " + str(res) + ", FSI tolerance = " + str(tol)
    #    raise Exception(ccolors.ANSI_RED + "FSI algo failed to converge!" + ccolors.ANSI_RESET)
    tests = CTests()
    tests.add(CTest('Lift coefficient', resultA[2], -0.425, 1e-1, False)) # rel. tol. of 10%
    tests.add(CTest('Drag coefficient', resultA[3], 3.124, 1e-1, False)) # rel. tol. of 10%
    tests.add(CTest('Displacement (180, TZ)', resultS[3], -0.000581, 1e-1, False)) # rel. tol. of 10%
    tests.add(CTest('Displacement (181, TZ)', resultS[2], -0.000810, 1e-1, False)) # rel. tol. of 10%
    tests.run()

def getFsiP():
    """Fsi parameters"""
    import os
    filePath = os.path.abspath(os.path.dirname(__file__))
    p = {}
    # Solvers and config files
    p['fluidSolver'] = 'SU2'
    p['solidSolver'] = 'SU2'
    p['cfdFile'] = os.path.join(filePath, 'config_channel.cfg')
    p['csdFile'] = os.path.join(filePath, 'config_cantilever.cfg')
    p['computation'] = 'Direct'
    # FSI objects
    p['interpolator'] = 'Matching'
    p['criterion'] = 'Displacements'
    p['algorithm'] = 'StaticBGS'
    # FSI parameters
    p['compType'] = 'steady'
    p['nDim'] = 2
    p['dt'] = 0.
    p['tTot'] = 0.05
    p['timeItTresh'] = -1
    p['tol'] = 1e-8
    p['maxIt'] = 15
    p['omega'] = 1.0
    p['rbfRadius'] = .3
    p['nodalLoadsType'] = 'force'
    return p

def main():
    import cupydo.interfaces.Cupydo as cupy
    p = getFsiP() # get parameters
    cupydo = cupy.CUPyDO(p) # create fsi driver
    cupydo.run() # run fsi process
    test(cupydo.algorithm.errValue, p['tol']) # check the results
    # eof
    #print('')

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':
    main()