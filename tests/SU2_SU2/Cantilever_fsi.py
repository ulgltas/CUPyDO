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
    with open("SolidSolution.ascii", 'rb') as f:
        lines = f.readlines()
    resultS = np.genfromtxt(lines[-1:], delimiter=None)

    # Check convergence and results
    # residual for test is 0.002
    #if (res > tol):
    #    print "\n\n" + "FSI residual = " + str(res) + ", FSI tolerance = " + str(tol)
    #    raise Exception("FSI algo failed to converge!")
    tests = CTests()
    tests.add(CTest('Lift coefficient', resultA[2], -0.425, 1e-1, False)) # rel. tol. of 10%
    tests.add(CTest('Drag coefficient', resultA[3], 3.124, 1e-1, False)) # rel. tol. of 10%
    tests.add(CTest('Displacement (110, Y)', resultS[6], -0.000581, 1e-1, False)) # rel. tol. of 10%
    tests.add(CTest('Displacement (100, Y)', resultS[3], -0.000810, 1e-1, False)) # rel. tol. of 10%
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
    p['computation'] = 'direct'
    # FSI objects
    p['interpolator'] = 'matching'
    p['interpType'] = 'conservative'
    p['algorithm'] = 'staticBGS'
    # FSI parameters
    p['regime'] = 'steady'
    p['nDim'] = 2
    p['dt'] = 0.
    p['tTot'] = 0.05
    
    p['dtSave'] = 0
    p['tol'] = 1e-8
    p['maxIt'] = 15
    p['omega'] = 1.0
    p['rbfRadius'] = .3
    # Solid parameters
    p['extractors'] = [2, 3]
    p['surfaceFilename'] = 'surface_solid'
    p['surfaceExtension'] = 'vtu'
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
    #print('')

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':
    main()
