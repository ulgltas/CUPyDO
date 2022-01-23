import cupydo.interfaces.Cupydo as cupy
import os

# %% Input Parameters

def getFsiP():

    p = dict()
    path = os.path.abspath(os.path.dirname(__file__))
    p['cfdFile'] = path+'\input_pfem.lua'

    # Metafor and PFEM solvers
    
    p['fluidSolver'] = 'Pfem3D'
    p['solidSolver'] = 'Metafor'
    p['csdFile'] = 'input_meta'
    
    # FSI objects

    p['criterion'] = 'Displacements'
    p['interpolator'] = 'Matching'
    p['algorithm'] = 'IQN_ILS'
    
    # FSI parameters

    p['firstItTgtMat'] = False
    p['computation'] = 'direct'
    p['compType'] = 'unsteady'
    p['timeItTresh'] = 0
    p['dtSave'] = 0.001
    p['omega'] = 0.5
    p['nSteps'] = 20
    p['maxIt'] = 25
    p['tTot'] = 10
    p['tol'] = 1e-6
    p['dt'] = 1e-3
    p['nDim'] = 2
    
    return p

# %% Main Function

def main():

    param = getFsiP()
    cupydo = cupy.CUPyDO(param)
    cupydo.run()

if __name__=='__main__':
    main()