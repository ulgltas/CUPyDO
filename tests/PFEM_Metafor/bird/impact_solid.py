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


from wrap import *

metafor = None


def params(q={}):
    """ default model parameters
    """
    p = {}
    p['tolNR'] = 1.0e-4        # Newton-Raphson tolerance
    p['tend'] = 2.            # final time
    p['dtmax'] = 0.005         # max time step
    p['bndno'] = 13            # interface boundary number
    p['saveAllFacs'] = False # keep the Fac corresponding to the end of the time step

    # BC type
    # p['bctype']     = 'pressure'     # uniform pressure
    # p['bctype']     = 'deadload'     # uniform nodal load
    # p['bctype']     = 'pydeadload1'  # uniform nodal load (python)
    p['bctype'] = 'pydeadloads'  # variable loads
    # p['bctype']     = 'slave'     # variable loads (mpi)

    p['exporter'] = Extractor()

    p.update(q)
    return p


def getMetafor(p={}):
    global metafor
    if metafor:
        return metafor
    metafor = Metafor()

    p = params(p)

    domain = metafor.getDomain()
    geometry = domain.getGeometry()
    geometry.setDimPlaneStrain(1.0)

    # import .geo
    from toolbox.gmshOld import GmshImport
    f = os.path.join(os.path.dirname(__file__), "impact.msh")
    importer = GmshImport(f, domain)
    importer.execute2D()

    groupset = domain.getGeometry().getGroupSet()

    # solid elements / material
    interactionset = domain.getInteractionSet()

    app1 = FieldApplicator(1)
    app1.push(groupset(17))  # physical group 17: panneau
    interactionset.add(app1)

    materset = domain.getMaterialSet()
    materset.define(1, EvpIsoHHypoMaterial)
    mater1 = materset(1)
    mater1.put(MASS_DENSITY,    2700.0)  # [kg/m³]
    mater1.put(ELASTIC_MODULUS, 69.0e9)  # [Pa]
    mater1.put(POISSON_RATIO,   0.3)   # [-]
    mater1.put(YIELD_NUM, 1)

    materlawset = domain.getMaterialLawSet()
    materlawset.define(1, LinearIsotropicHardening)
    materlawset(1).put(IH_SIGEL, 300.0e6)  # [Pa]
    materlawset(1).put(IH_H,    1000.0e6)  # [Pa]

    prp = ElementProperties(Volume2DElement)
    app1.addProperty(prp)
    prp.put(MATERIAL, 1)
    prp.put(CAUCHYMECHVOLINTMETH, VES_CMVIM_SRIPR)

    # boundary conditions
    loadingset = domain.getLoadingSet()

    # Physical Line(101) - clamped side of the beam
    x_ext_left = -0.4
    y_ext_left = -0.00635
    x_ext_right = 0.4
    y_ext_right = y_ext_left
    delta = 0.00001
    groupset.add(Group(50))
    groupset(50).addMeshPointsInBox(x_ext_left-delta, x_ext_left +
                                    delta, y_ext_left-delta, y_ext_left+delta, -delta, delta)
    groupset(50).addMeshPointsInBox(x_ext_right-delta, x_ext_right +
                                    delta, y_ext_right-delta, y_ext_right+delta, -delta, delta)
    loadingset.define(groupset(50), Field1D(TY, RE), 0.)
    loadingset.define(groupset(14), Field1D(TX, RE))
    # Physical Line(102) - free surface of the beam
    # Physical Line(103) - upper surface of the beam (for tests only)

    mim = metafor.getMechanicalIterationManager()
    mim.setResidualTolerance(p['tolNR'])
    # mim.setResidualComputationMethod(Method4ResidualComputation(1000.))

    ti = AlphaGeneralizedTimeIntegration(metafor)
    metafor.setTimeIntegration(ti)

    # visu
    if 0:
        tsm = metafor.getTimeStepManager()
        tsm.setInitialTime(0.0, 1.0)
        tsm.setNextTime(1.0, 1, 1.0)

    # results (DOES NOT WORK PROPERLY!)
    vmgr = metafor.getValuesManager()
    # vmgr.add(1, TdFieldValueExtractor(metafor, groupset(17), THERMODYN_TRAV_FINT), 'panel_work_int_F')
    vmgr.add(2, TdFieldValueExtractor(
        metafor, groupset(17), THERMODYN_EN_CIN), 'panel_kin_En')
    vmgr.add(3, TdFieldValueExtractor(
        metafor, groupset(17), THERMODYN_EN_DIS), 'panel_diss_En')
    vmgr.add(4, TdFieldValueExtractor(metafor, groupset(
        17), THERMODYN_TRAV_FEXT), 'panel_work_ext_F')

    metafor.getDbTdFieldsValueManager().setComputeEXW(True)
    metafor.getDbTdFieldsValueManager().setComputeINW(True)

    return metafor

class Extractor(object):
    def __init__(self):

        self.metafor = metafor
        

    def write(self,extractor):

        file = open(extractor.buildName()+'.ascii', 'a')
        
        file.write('{0:12.6f}\t'.format(self.metafor.getCurrentTime()))
        file.write('{0:12.6f}\n'.format(extractor.extract()[0]))
        file.close()

    def execute(self):

        groupset = self.metafor.getDomain().getGeometry().getGroupSet()
        self.write(TdFieldValueExtractor(metafor, groupset(17), THERMODYN_TRAV_FINT))
        self.write(DbNodalValueExtractor(groupset(17), Field1D(TY, RE), SortByDist0(0,0,0), 1))
        