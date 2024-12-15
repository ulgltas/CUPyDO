import cupydo.interfaces.Cupydo as cupy
import os

# Input Parameters

def getFsiP():

    p = dict()
    path = os.path.abspath(os.path.dirname(__file__))
    p['cfdFile'] = path+'/input_pfem.lua'

    # Metafor and PFEM solvers
    
    p['fluidSolver'] = 'Pfem3D'
    p['solidSolver'] = 'Metafor'
    p['csdFile'] = 'input_meta'
    
    # FSI objects

    p['interpolator'] = 'matching'
    p['interpType'] = 'consistent'
    p['algorithm'] = 'IQN_MVJ'
    
    # FSI parameters

    p['firstItTgtMat'] = False
    p['computation'] = 'direct'
    p['regime'] = 'unsteady'
    
    p['omega'] = 0.5
    p['dtSave'] = 1e-3
    p['maxIt'] = 25
    p['tTot'] = 2
    p['dt'] = 1e-4
    p['criterion'] = 'relative'
    p['nDim'] = 3
    p['qrFilter'] = None
    
    # Coupling Type

    p['mechanical'] = True
    p['mechanicalTol'] = 1e-8
    p['thermal'] = False
    return p

# Main Function

def main():

    param = getFsiP()
    cupydo = cupy.CUPyDO(param)
    cupydo.run()

    # Compare the results with a reference solution

    base = os.path.dirname(__file__)
    os.system(f'python3 {base}/battery.py')

    # Count the number of iterations

    n_iter = cupydo.algorithm.totNbOfFSIIt
    n_step = cupydo.algorithm.step.timeIter

    print('Iterations =', n_iter)
    print('Time steps =', n_step)
    print('Mean =', n_iter/n_step)

if __name__=='__main__':
    main()