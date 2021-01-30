#!/usr/bin/env python
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
    from cupydo.testing import *
    # Read results from file
    with open("AerodynamicCoeff.ascii", 'rb') as f:
        lines = f.readlines()
    resultA = np.genfromtxt(lines[-1:], delimiter=None)

    # Check convergence and results
    if (res > tol):
        print "\n\n" + "FSI residual = " + str(res) + ", FSI tolerance = " + str(tol)
        raise Exception(ccolors.ANSI_RED + "FSI algo failed to converge!" + ccolors.ANSI_RESET)
    tests = CTests()
    tests.add(CTest('Lift coefficient', resultA[2], -0.000010, 1e-1, False)) # rel. tol. of 10%
    tests.add(CTest('Drag coefficient', resultA[3], 1.7198, 1e-1, False)) # rel. tol. of 10%
    tests.run()

def getFsiP():
    import os
    filePath = os.path.abspath(os.path.dirname(__file__))
    p = {}
    # Solvers and config files
    p['fluidSolver'] = 'SU2'
    p['solidSolver'] = 'GetDP'
    p['cfdFile'] = os.path.join(filePath, 'CHT_Cylinder_FFTB_SU2Conf.cfg')
    p['csdFile'] = 'CHT_Cylinder_FFTB'
    # FSI objects
    p['interpolator'] = 'Matching'
    p['criterion'] = 'Displacements'
    p['algorithm'] = 'StaticBGS'
    # FSI parameters
    p['compType'] = 'steady'
    p['nDim'] = 2
    p['dt'] = 0.0
    p['tTot'] = 0.0
    p['timeItTresh'] = -1
    p['tol'] = 1e-6
    p['maxIt'] = 1000
    p['omega'] = [1.0, 1.0]
    p['nodalLoadsType'] = 'pressure'
    #Param for CHT and GETDP
    p['tolCHT'] = 0.1
    p['chtScheme'] = 'FFTB'
    p['getDPExtractNode'] = 7
    p['GetDPResolution'] = 'HeatEquation'
    p['mechanicalCoupling'] = False
    p['thermalcoupling'] = True
    p['fluidTWallInit'] = 330.0
    return p

def main():

    import cupydo.interfaces.Cupydo as cupy
    p = getFsiP()
    cupydo = cupy.CUPyDO(p)
    cupydo.run()
    test(cupydo.algorithm.errValue_CHT, p['tolCHT']) # check the results

    # eof
    print ''

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':
    main()
