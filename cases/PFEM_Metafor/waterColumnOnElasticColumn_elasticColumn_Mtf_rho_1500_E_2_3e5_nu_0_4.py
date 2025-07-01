# -*- coding: latin-1; -*-
# bends a simple beam from a gmsh file

from wrap import *

metafor = None

def getMetafor(p={}):
    global metafor
    if metafor: return metafor
    metafor = Metafor()

    domain = metafor.getDomain()
    geometry = domain.getGeometry()
    geometry.setDimPlaneStrain(1.0)

    # import .geo
    from toolbox.gmshOld import GmshImport
    f = os.path.join(os.path.dirname(__file__), "waterColumnOnElasticColumn_Mtf_Pfem.msh")
    importer = GmshImport(f, domain)
    importer.execute2D()

    groupset = domain.getGeometry().getGroupSet()
    p['FSI'] = groupset(14)

    # solid elements / material
    interactionset = domain.getInteractionSet()

    app1 = FieldApplicator(1)
    app1.push( groupset(17) )  # physical group 100: beam
    interactionset.add( app1 )

    materset = domain.getMaterialSet()
    materset.define( 1, ElastHypoMaterial )
    mater1 = materset(1)
    mater1.put(MASS_DENSITY,    1500.0)  # [kg/m³]
    mater1.put(ELASTIC_MODULUS, 2.3e5)  # [Pa]
    mater1.put(POISSON_RATIO,   0.4)   # [-]

    prp = ElementProperties(TriangleVolume2DElement)
    app1.addProperty(prp)
    prp.put (MATERIAL, 1)
    prp.put(GRAVITY_Y, -10) # m/s2
    # prp.put(CAUCHYMECHVOLINTMETH,VES_CMVIM_SRIPR)
    
    # boundary conditions
    loadingset = domain.getLoadingSet()

    #Physical Line(101) - clamped side of the beam
    loadingset.define(groupset(12), Field1D(TX,RE))
    loadingset.define(groupset(12), Field1D(TY,RE))
    loadingset.define(groupset(13), Field1D(TX,RE))
    #Physical Line(102) - free surface of the beam  
    #Physical Line(103) - upper surface of the beam (for tests only)


    mim = metafor.getMechanicalIterationManager()
    mim.setMaxNbOfIterations(4)
    mim.setResidualTolerance(1.0e-6)

    ti = AlphaGeneralizedTimeIntegration(metafor)
    ti.setAlphaM(0.)
    ti.setAlphaF(0.)
    metafor.setTimeIntegration(ti)

    # visu
    if 0:
        tsm = metafor.getTimeStepManager()
        tsm.setInitialTime(0.0, 1.0)
        tsm.setNextTime(1.0, 1, 1.0)

    # results
    vmgr = metafor.getValuesManager()
    vmgr.add(1, MiscValueExtractor(metafor, EXT_T), 'time_mtf')
    vmgr.add(2, DbNodalValueExtractor(groupset(14), Field1D(TY,RE)), 'dy')
    
    p['exporter'] = None
    return metafor




