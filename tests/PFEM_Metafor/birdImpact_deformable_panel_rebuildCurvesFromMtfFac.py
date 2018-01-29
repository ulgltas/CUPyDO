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

from wrap import *
import math
import os, os.path, sys, time, string
import toolbox

module = toolbox.utilities.load('birdImpact_deformable_panel_panel_steel_Mtf')
toolbox.utilities.setTheWDir("D:/fsi/tests/workspace/birdImpact_deformable_panel_Mtf_Pfem")

metafor = module.getMetafor()

valuesmanager = metafor.getValuesManager()
valuesmanager.clear()
valuesmanager.add(1, TdFieldValueExtractor(metafor, metafor.getDomain().getGeometry().getGroupSet()(17), THERMODYN_EN_CIN), 'panel_kin_En')
valuesmanager.add(2, TdFieldValueExtractor(metafor, metafor.getDomain().getGeometry().getGroupSet()(17), THERMODYN_EN_DIS), 'panel_diss_En')
valuesmanager.add(3, DbNodalValueExtractor(metafor.getDomain().getGeometry().getGroupSet()(17), Field1D(INW, AB)), SumOperator(), 'panel_work_int_F')
valuesmanager.add(4, DbNodalValueExtractor(metafor.getDomain().getGeometry().getGroupSet()(17), Field1D(EXW, AB)), SumOperator(), 'panel_work_ext_F') 

toolbox.utilities.rebuildCurves()