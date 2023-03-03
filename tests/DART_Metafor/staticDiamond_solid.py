#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from wrap import *
from wrap.mtCompositesw import *
from math import *

def params(_p):
    p={}
    p['relTol'] = 1e-6
    p['maxIt'] = 20
    p['bndno'] = 111
    p['saveAllFacs'] = False
                       
    p.update(_p)
    return p

def getMetafor(p={}):
    metafor = Metafor()
    p = params(p)

    # -- Geometry container and computation type
    domain = metafor.getDomain()
    geometry = domain.getGeometry()
    geometry.setDimPlaneStrain(1.0)

    # -- Import the geometry
    from toolbox.gmshOld import GmshImport
    f = os.path.join(os.path.dirname(__file__), "models/diamond_solid.geo")
    importer = GmshImport(f, domain)
    importer.execute2D()    
    groupset = domain.getGeometry().getGroupSet()

    # -- Define the material
    matset = domain.getMaterialSet()
    matset.define( 1, ElastHypoMaterial )
    mat1 = matset(1)
    mat1.put(MASS_DENSITY,    100.0)  # [kg/m³]
    mat1.put(ELASTIC_MODULUS, 1e8)  # [Pa]
    mat1.put(POISSON_RATIO,   0.35)   # [-]

    # -- Define quad solid elements
    interactionset = domain.getInteractionSet()
    fieldapp1 = FieldApplicator(1)
    fieldapp1.push( groupset(121) )    #physical group 'material' = volumic mesh of the wing
    interactionset.add( fieldapp1 )
    prp = ElementProperties(Volume2DElement)
    fieldapp1.addProperty(prp)
    prp.put(MATERIAL, 1)
    prp.put(CAUCHYMECHVOLINTMETH,VES_CMVIM_STD)

    # -- Define tri solid elements
    fieldapp2 = FieldApplicator(2)
    fieldapp2.push( groupset(121) )
    interactionset.add( fieldapp2 )
    prp2 = ElementProperties(TriangleVolume2DElement)
    fieldapp2.addProperty(prp2)
    prp2.put(MATERIAL, 1)
    prp2.put(CAUCHYMECHVOLINTMETH,VES_CMVIM_STD)

    # -- Boundary conditions
    loadingset = domain.getLoadingSet()
    # Clamped node
    loadingset.define(groupset(103), Field1D(TX,RE))
    loadingset.define(groupset(103), Field1D(TY,RE))
    # Supported node
    loadingset.define(groupset(102), Field1D(TX,RE))

    # -- Numerical parameters
    mim = metafor.getMechanicalIterationManager()
    mim.setMaxNbOfIterations(p['maxIt'])
    mim.setResidualTolerance(p['relTol'])

    # -- for pure Metafor testing
    if 0:         
        # -- Time integration
        ti = QuasiStaticTimeIntegration(metafor)
        metafor.setTimeIntegration(ti)

        # -- Value manager
        vmgr = metafor.getValuesManager()
        vmgr.add(1, DbNodalValueExtractor(groupset(101), Field1D(TY,RE)), 'TE_dy')
        
        # prescribe node displacement
        fct = PieceWiseLinearFunction()
        fct.setData(0.0, 0.0)
        fct.setData(1.0, 1.0)
        loadingset.define(groupset(101), Field1D(TY,RE), 0.2, fct)

        tsm = metafor.getTimeStepManager()
        tsm.setInitialTime(0.0, 0.1)
        tsm.setNextTime(1.0, 1, 0.1)
  
    return metafor

def getRealTimeExtractorsList(Mtf):

    extractorsList = list()
    groupset = Mtf.getDomain().getGeometry().getGroupSet()

    # --- Extractors list starts --- #
    extractor = DbNodalValueExtractor(groupset(101), Field1D(TY,RE))
    extractorsList.append(extractor)
    # --- Extractors list ends --- #

    return extractorsList