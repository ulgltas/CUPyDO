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
    f = os.path.join(os.path.dirname(__file__), "MovingSquareCoarse.msh")
    importer = GmshImport(f, domain)
    importer.execute2D()

    groupset = domain.getGeometry().getGroupSet()
    p['FSI'] = groupset(3)

    # solid elements / material
    interactionset = domain.getInteractionSet()

    app1 = FieldApplicator(1)
    app1.push( groupset(5) )  # physical group 5: square
    interactionset.add( app1 )
    
    materset = domain.getMaterialSet()
    materset.define( 1, ElastHypoMaterial )
    mater1 = materset(1)
    mater1.put(MASS_DENSITY,    1.2)  # [kg/m3]
    mater1.put(ELASTIC_MODULUS, 1e9)  # [Pa] (Value given by Franci et al. (CMAME, 2016) --> to be verified)
    mater1.put(POISSON_RATIO,   0.49)   # [-]

    prp = ElementProperties(Volume2DElement)
    app1.addProperty(prp)
    prp.put (MATERIAL, 1)
    prp.put(GRAVITY_Y, -10.) # m/s2
    prp.put(CAUCHYMECHVOLINTMETH,VES_CMVIM_STD)
    
    mim = metafor.getMechanicalIterationManager()
    mim.setMaxNbOfIterations(10)
    mim.setResidualTolerance(1.0e-6)

    ti = AlphaGeneralizedTimeIntegration(metafor)
    '''ti.setAlphaM(0.)
    ti.setAlphaF(0.)
    ti.setBeta0(0.25)
    ti.setGamma0(0.5)'''
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




