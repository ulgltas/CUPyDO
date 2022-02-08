#!/usr/bin/env python
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
    geometry.setDim3D()

    # -- Import the geometry
    from toolbox.gmsh import GmshImport
    f = os.path.join(os.path.dirname(__file__), "models/agard445_solid.geo")
    importer = GmshImport(f, domain)
    importer.execute()    
    groupset = domain.getGeometry().getGroupSet()

    # -- Define the material
    matset = domain.getMaterialSet()
    matset.define( 1, ElastOrthoHypoMaterial )
    matset.define( 2, ElastOrthoHypoMaterial )
    mat = matset(1)
    mat.put(MASS_DENSITY, 381.98)          # [kg/m3]
    mat.put(YOUNG_MODULUS_1, 3.151e9)      # [Pa]
    mat.put(YOUNG_MODULUS_2, 0.4162e9)     # [Pa]
    mat.put(YOUNG_MODULUS_3, 0.4162e9)     # [Pa]
    mat.put(POISSON_RATIO_12, 0.31)        # [-]
    mat.put(POISSON_RATIO_13, 0.31)        # [-]
    mat.put(POISSON_RATIO_23, 0.31)        # [-]
    mat.put(SHEAR_MODULUS_12, 0.4392e9)    # [Pa]
    mat.put(SHEAR_MODULUS_13, 0.4392e9)    # [Pa]
    mat.put(SHEAR_MODULUS_23, 0.4392e9)    # [Pa]
    sweepOrtho = pi/4.0
    mat.put(ORTHO_AX1_X, sin(sweepOrtho))
    mat.put(ORTHO_AX1_Y, cos(sweepOrtho))
    mat.put(ORTHO_AX1_Z, 0.0)
    mat.put(ORTHO_AX2_X, cos(sweepOrtho))
    mat.put(ORTHO_AX2_Y, -sin(sweepOrtho))
    mat.put(ORTHO_AX2_Z, 0.0)

    mat2 = matset(2)
    mat2.put(MASS_DENSITY, 381.98)          # [kg/m3]
    mat2.put(YOUNG_MODULUS_1, 3.151e9)      # [Pa]
    mat2.put(YOUNG_MODULUS_2, 0.4162e9)     # [Pa]
    mat2.put(YOUNG_MODULUS_3, 0.4162e9)     # [Pa]
    mat2.put(POISSON_RATIO_12, 0.31)        # [-]
    mat2.put(POISSON_RATIO_13, 0.31)        # [-]
    mat2.put(POISSON_RATIO_23, 0.31)        # [-]
    mat2.put(SHEAR_MODULUS_12, 0.4392e9)    # [Pa]
    mat2.put(SHEAR_MODULUS_13, 0.4392e9)    # [Pa]
    mat2.put(SHEAR_MODULUS_23, 0.4392e9)    # [Pa]
    sweepOrtho = -pi/4.0
    mat2.put(ORTHO_AX1_X, sin(sweepOrtho))
    mat2.put(ORTHO_AX1_Y, cos(sweepOrtho))
    mat2.put(ORTHO_AX1_Z, 0.0)
    mat2.put(ORTHO_AX2_X, cos(sweepOrtho))
    mat2.put(ORTHO_AX2_Y, -sin(sweepOrtho))
    mat2.put(ORTHO_AX2_Z, 0.0)

    # -- Define hexa solid elements
    interactionset = domain.getInteractionSet()
    fieldapp1 = FieldApplicator(1)
    fieldapp1.push( groupset(101) )    #physical group 'material' = volumic mesh of the wing
    interactionset.add( fieldapp1 )
    prp = ElementProperties(Volume3DElement)
    fieldapp1.addProperty(prp)
    prp.put(MATERIAL, 1)
    prp.put(CAUCHYMECHVOLINTMETH,VES_CMVIM_EAS)
    prp.put(PEAS, 1e-10)

    # -- Define prism solid elements
    fieldapp2 = FieldApplicator(2)
    fieldapp2.push( groupset(101) )
    interactionset.add( fieldapp2 )
    prp2 = ElementProperties(PentaVolume3DElement)
    fieldapp2.addProperty(prp2)
    prp2.put(MATERIAL, 1)
    prp2.put(CAUCHYMECHVOLINTMETH,VES_CMVIM_STD)

    # -- Define hexa solid elements: symmetry
    interactionset = domain.getInteractionSet()
    fieldapp3 = FieldApplicator(3)
    fieldapp3.push( groupset(102) )    #physical group 'material' = volumic mesh of the wing
    interactionset.add( fieldapp3 )
    prp3 = ElementProperties(Volume3DElement)
    fieldapp3.addProperty(prp3)
    prp3.put(MATERIAL, 2)
    prp3.put(CAUCHYMECHVOLINTMETH,VES_CMVIM_EAS)
    prp3.put(PEAS, 1e-10)

    # -- Define prism solid elements: symmetry
    fieldapp4 = FieldApplicator(4)
    fieldapp4.push( groupset(102) )
    interactionset.add( fieldapp4 )
    prp4 = ElementProperties(PentaVolume3DElement)
    fieldapp4.addProperty(prp4)
    prp4.put(MATERIAL, 2)
    prp4.put(CAUCHYMECHVOLINTMETH,VES_CMVIM_STD)

    # -- Boundary conditions
    loadingset = domain.getLoadingSet()
    # Clamped face
    loadingset.define(groupset(112), Field1D(TX,RE))
    loadingset.define(groupset(112), Field1D(TY,RE))
    loadingset.define(groupset(112), Field1D(TZ,RE))

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
        vmgr.add(1, DbNodalValueExtractor(groupset(121), Field1D(TZ,RE)), 'TE_dz')
        vmgr.add(2, DbNodalValueExtractor(groupset(122), Field1D(TZ,RE)), 'LE_dz')
        
        # prescribe node displacement
        fct = PieceWiseLinearFunction()
        fct.setData(0.0, 0.0)
        fct.setData(1.0, 1.0)
        loadingset.define(groupset(122), Field1D(TZ,RE), 0.2, fct)

        tsm = metafor.getTimeStepManager()
        tsm.setInitialTime(0.0, 0.1)
        tsm.setNextTime(1.0, 1, 0.1)
  
    return metafor

def getRealTimeExtractorsList(Mtf):

    extractorsList = list()
    groupset = Mtf.getDomain().getGeometry().getGroupSet()

    # --- Extractors list starts --- #
    extractor0 = DbNodalValueExtractor(groupset(121), Field1D(TZ,RE))
    extractorsList.append(extractor0)
    extractor1 = DbNodalValueExtractor(groupset(122), Field1D(TZ,RE))
    extractorsList.append(extractor1)
    # --- Extractors list ends --- #

    return extractorsList