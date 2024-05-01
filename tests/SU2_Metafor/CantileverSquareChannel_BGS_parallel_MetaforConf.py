#! /usr/bin/env python3
# -*- coding: utf8 -*-

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
    p={}
    p['tolNR']      = 1.0e-7        # Newton-Raphson tolerance
    p['bndno'] = 102
    p['exporter'] = Extractor()
                                       
    p.update(q)
    return p

def getMetafor(p={}):
    global metafor
    if metafor: return metafor
    metafor = Metafor()

    p = params(p)

    domain = metafor.getDomain()
    geometry = domain.getGeometry()
    geometry.setDimPlaneStrain(1.0)

    # import .geo
    from toolbox.gmshOld import GmshImport
    f = os.path.join(os.path.dirname(__file__), "models/CantileverSquareChannel_BGS_parallel_solidMesh.geo")
    importer = GmshImport(f, domain)
    importer.execute2D()

    groupset = domain.getGeometry().getGroupSet()    

    # solid elements / material
    interactionset = domain.getInteractionSet()

    app1 = FieldApplicator(1)
    app1.push( groupset(100) )  # physical group 100 = meshed beam
    interactionset.add( app1 )

    materset = domain.getMaterialSet()
    materset.define( 1, ElastHypoMaterial )
    mater1 = materset(1)
    mater1.put(MASS_DENSITY,    100.0)  # [kg/m³]
    mater1.put(ELASTIC_MODULUS, 2.5e5)  # [Pa]
    mater1.put(POISSON_RATIO,   0.35)   # [-]

    prp = ElementProperties(Volume2DElement)
    app1.addProperty(prp)
    prp.put (MATERIAL, 1)
    #prp.put(CAUCHYMECHVOLINTMETH,VES_CMVIM_STD)
    prp.put(CAUCHYMECHVOLINTMETH,VES_CMVIM_SRI)
    
    # boundary conditions
    loadingset = domain.getLoadingSet()

    #Physical Line(101) - clamped side of the beam
    loadingset.define(groupset(101), Field1D(TX,RE))
    loadingset.define(groupset(101), Field1D(TY,RE))


    mim = metafor.getMechanicalIterationManager()
    mim.setMaxNbOfIterations(10)
    mim.setResidualTolerance(p['tolNR'])

    ti = AlphaGeneralizedTimeIntegration(metafor)
    metafor.setTimeIntegration(ti)

    return metafor

class Extractor(object):
    def __init__(self):

        self.metafor = metafor
        

    def to_ascii(self,extractor):

        file = open(extractor.buildName()+'.ascii', 'a')
        
        file.write('{0:12.6f}\t'.format(self.metafor.getCurrentTime()))
        file.write('{0:12.6f}\n'.format(extractor.extract()[0]))
        file.close()

    def write(self):

        groupset = self.metafor.getDomain().getGeometry().getGroupSet()
        self.to_ascii(DbNodalValueExtractor(groupset(104), Field1D(TY,RE)))
        


