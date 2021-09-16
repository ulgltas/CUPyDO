import cupydo.interfaces.Cupydo as cupy

# %% Input Parameters

def getFsiP():

    param = {}
    
    # Metafor and PFEM solvers
    
    param['fluidSolver'] = 'Pfem'
    param['solidSolver'] = 'Metafor'
    param['cfdFile'] = 'input_pfem'
    param['csdFile'] = 'input_meta'
    
    # FSI objects

    param['interpolator'] = 'Matching'
    param['criterion'] = 'Displacements'
    param['algorithm'] = 'AitkenBGS'
    
    # FSI parameters

    param['compType'] = 'unsteady'
    param['computation'] = 'direct'
    param['nDim'] = 2
    param['dt'] = 0.001
    param['tTot'] = 0.3
    param['timeItTresh'] = 0
    param['tol'] = 1e-6
    param['maxIt'] = 25
    param['omega'] = 0.5
    
    return param

param = getFsiP()
cupydo = cupy.CUPyDO(param)
cupydo.run()