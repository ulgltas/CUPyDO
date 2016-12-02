import os, sys

filePath = os.path.abspath(os.path.dirname(sys.argv[0]))
fileName = os.path.splitext(os.path.basename(__file__))[0]
fsiPath = '..'
sys.path.append(fsiPath)

import fsi
import MtfInterface

def main():
    
    fsi.load(fileName)
    
    U0 = 100
    N = 10
    R = 0.01
    d = 2.5*(R/N)
    
    t1 = 0.0
    dt = 2e-6
    tTot = 40*((4*R)/U0 + d/U0)
    
    solid = MtfInterface.MtfSolver('birdImpact_deformable_panel_panel_steel_Mtf')
    fluid = fsi.PfemSolver('birdImpact_deformable_panel_bird_Pfem', 13, dt)
    
    toll = 1.0e-6
    fsi_criterion = fsi.DispResidualBasedCriterion(toll)
    
    omega_max = 0.5
    
    fsi_algo = fsi.FixedPointAitkenRelaxationAlgorithm(solid, fluid, dt, tTot, fsi_criterion, omega_max)
    fsi_algo.run()
    
if __name__ == "__main__":
    main()