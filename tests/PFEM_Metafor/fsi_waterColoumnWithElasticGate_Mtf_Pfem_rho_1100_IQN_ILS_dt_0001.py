import os, sys

filePath = os.path.abspath(os.path.dirname(sys.argv[0]))
fileName = os.path.splitext(os.path.basename(__file__))[0]
FSICouplerPath = os.path.join(os.path.dirname(__file__), '../../')
appsPath = os.path.join(os.path.dirname(__file__), '../../apps')

sys.path.append(FSICouplerPath)
sys.path.append(appsPath)

from math import *
from optparse import OptionParser
import FSICoupler

def main(nogui): # NB, the argument 'nogui' is specific to PFEM only!
    
    p = {}
    p['tTot'] = 0.05
    p['mtfSaveAllFacs'] = False
    p['saveFreqPFEM'] = 1000
    
    from PFEM_Metafor.fsi_waterColoumnWithElasticGate_Mtf_Pfem_rho_1100_IQN_ILS_dt_0001 import main
    main(p, nogui)

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':
    
    parser=OptionParser()
    parser.add_option("--nogui", action="store_true",
                        help="Specify if we need to use the GUI", dest="nogui", default=False)

    (options, args)=parser.parse_args()
    
    nogui = options.nogui
    
    main(nogui)