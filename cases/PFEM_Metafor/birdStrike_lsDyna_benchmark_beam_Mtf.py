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
    f = os.path.join(os.path.dirname(__file__), "birdStrike_lsDyna_benchmark_Mtf_Pfem.msh")
    importer = GmshImport(f, domain)
    importer.execute2D()

    groupset = domain.getGeometry().getGroupSet()
    p['FSI'] = groupset(15)  

    # solid elements / material
    interactionset = domain.getInteractionSet()

    app1 = FieldApplicator(1)
    app1.push( groupset(16) )  # physical group 100: beam
    interactionset.add( app1 )

    materset = domain.getMaterialSet()
    materset.define( 1, ElastHypoMaterial )
    mater1 = materset(1)
    mater1.put(MASS_DENSITY,    1000.0)  # [kg/m³]
    mater1.put(ELASTIC_MODULUS, 1.0e3)  # [Pa]
    mater1.put(POISSON_RATIO,   0.)   # [-] 

    prp = ElementProperties(Volume2DElement)
    app1.addProperty(prp)
    prp.put (MATERIAL, 1)
    prp.put(CAUCHYMECHVOLINTMETH,VES_CMVIM_STD)
    '''prp.put(EASS, 2)
    prp.put(EASV, 2)'''
    
    # boundary conditions
    loadingset = domain.getLoadingSet()

    #Physical Line(101) - clamped side of the beam
    loadingset.define(groupset(13), Field1D(TX,RE))
    loadingset.define(groupset(13), Field1D(TY,RE))
    #Physical Line(102) - free surface of the beam  
    #Physical Line(103) - upper surface of the beam (for tests only)


    mim = metafor.getMechanicalIterationManager()
    mim.setResidualTolerance(1.0e-7)
    # mim.setResidualComputationMethod(Method4ResidualComputation(1000.)) 

    ti = AlphaGeneralizedTimeIntegration(metafor)
    metafor.setTimeIntegration(ti)

    # visu
    if 0:
        tsm = metafor.getTimeStepManager()
        tsm.setInitialTime(0.0, 1.0)
        tsm.setNextTime(1.0, 1, 1.0)

    # results
    #vmgr = metafor.getValuesManager()
    #vmgr.add(1, MiscValueExtractor(metafor, EXT_T), 'time')
    #vmgr.add(2, DbNodalValueExtractor(groupset(104), Field1D(TY,RE)), 'dy')
    
    p['exporter'] = None
    return metafor





