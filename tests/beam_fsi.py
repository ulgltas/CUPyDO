import os, sys

filePath = os.path.abspath(os.path.dirname(sys.argv[0]))
fileName = os.path.splitext(os.path.basename(__file__))[0]
fsiPath = '..'
sys.path.append(fsiPath)

import fsi

def main():
    
    fsi.load(fileName)
    
    t1 = 0.0
    
    solid = fsi.MtfSolver('beam', 103, t1)
    solid.saveAllFacs = False
    
    fsi_algo = fsi.fsiSolidTestAlgorithm(solid)
    fsi_algo.run()
    
if __name__ == "__main__":
    main()