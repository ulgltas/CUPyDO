#! /usr/bin/env python
# -*- coding: latin-1; -*-

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

import os, sys

filePath = os.path.abspath(os.path.dirname(sys.argv[0]))
fileName = os.path.splitext(os.path.basename(__file__))[0]
FSICouplerPath = os.path.join(os.path.dirname(__file__), '../../')
appsPath = os.path.join(os.path.dirname(__file__), '../../apps')

sys.path.append(FSICouplerPath)
sys.path.append(appsPath)

from math import *
from optparse import OptionParser

def main(nogui): # NB, the argument 'nogui' is specific to PFEM only!
    
    p = {}
    p['tTot'] = 1e-4
    p['mtfSaveAllFacs'] = False
    p['saveFreqPFEM'] = 1000
    
    from PFEM_Metafor.fsi_birdImpact_deformable_panel_alu_Mtf_Pfem_Smojver_Axisym import main
    main(p, nogui)

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':
    
    parser=OptionParser()
    parser.add_option("--nogui", action="store_true",
                        help="Specify if we need to use the GUI", dest="nogui", default=False)

    (options, args)=parser.parse_args()
    
    nogui = options.nogui
    
    main(nogui)