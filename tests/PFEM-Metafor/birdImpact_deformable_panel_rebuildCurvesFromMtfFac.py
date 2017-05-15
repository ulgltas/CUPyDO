# Reconstruct desired curves a posteriori from Metafor Fac files

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