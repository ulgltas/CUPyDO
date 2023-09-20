#! /usr/bin/env python
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


def test(res,tol,it):

    from cupydo.testing import CTest,CTests
    import numpy as np

    # Check convergence and results

    if (res > tol):
        print("\n\nFSI residual = "+str(res)+", FSI tolerance = "+str(tol))
        raise Exception("FSI algo failed to converge!")

    # Read results from file

    with open("Node_56_POS.ascii", 'rb') as f: lines = f.readlines()
    result = np.genfromtxt(lines[-1:], delimiter=None)

    tests = CTests()
    tests.add(CTest('Middle bar coordinate X', result[0], 0.5 , 1e-3, False))
    tests.add(CTest('Middle bar coordinate Y', result[1], -0.053139, 0.05, False))
    tests.add(CTest('Mean number of ISI iterations', it, 4, 1, True))
    tests.run()

# Input Parameters

def getFsiP():

    p = dict()

    # Metafor and PFEM solvers
    
    p['fluidSolver'] = 'Pfem'
    p['solidSolver'] = 'Metafor'
    p['cfdFile'] = 'beam_fluid'
    p['csdFile'] = 'beam_solid'
    
    # FSI objects

    p['interpolator'] = 'matching'
    p['interpType'] = 'conservative'
    p['algorithm'] = 'IQN_ILS'
    
    # FSI parameters

    p['firstItTgtMat'] = False
    p['computation'] = 'direct'
    p['regime'] = 'unsteady'
    
    p['dtSave'] = 0
    p['omega'] = 0.5
    p['maxIt'] = 25
    p['nSteps'] = 10
    p['tol'] = 1e-8
    p['dt'] = 0.1
    p['tTot'] = 20
    p['nDim'] = 2

    # Coupling Type

    p['mechanical'] = True
    p['thermal'] = False
    return p

# Main Function

def main():

    import cupydo.interfaces.Cupydo as cupy

    p = getFsiP() # get parameters
    cupydo = cupy.CUPyDO(p) # create fsi driver
    cupydo.run() # run fsi process

    cupydo.algorithm.FluidSolver.save(cupydo.algorithm.step.timeIter)
    test(cupydo.algorithm.errValue, p['tol'], cupydo.algorithm.getMeanNbOfFSIIt()) # check the results
    
    # eof
    print('')

if __name__=='__main__':
    main()