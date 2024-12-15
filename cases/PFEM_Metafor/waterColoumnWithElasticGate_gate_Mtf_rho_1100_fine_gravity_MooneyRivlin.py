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
    f = os.path.join(os.path.dirname(__file__), "waterColoumnWithElasticGate_Mtf_Pfem_fine.msh")
    importer = GmshImport(f, domain)
    importer.execute2D()

    groupset = domain.getGeometry().getGroupSet()
    p['FSI'] = groupset(17) 

    # solid elements / material
    interactionset = domain.getInteractionSet()

    app1 = FieldApplicator(1)
    app1.push( groupset(21) )  # physical group 100: beam
    interactionset.add( app1 )
    
    G = 2.4e6 #2.6e6
    mooneyrivlin_c1 = -1.2e6
    mooneyrivlin_c2 = G/2.0 - mooneyrivlin_c1
    penal = 1.3e6
    materset = domain.getMaterialSet()
    materset.define(1,MooneyRivlinHyperMaterial)
    mater1 = materset(1)
    mater1.put(MASS_DENSITY, 1100.)
    mater1.put(RUBBER_C1, mooneyrivlin_c1)
    mater1.put(RUBBER_C2, mooneyrivlin_c2)
    mater1.put(RUBBER_PENAL, penal)
    
    prp = ElementProperties(Volume2DElement)
    app1.addProperty(prp)
    prp.put (MATERIAL, 1)
    prp.put(GRAVITY_Y, -9.81) # m/s2
    prp.put(STIFFMETHOD, STIFF_NUMERIC)
    # prp.put(CAUCHYMECHVOLINTMETH,VES_CMVIM_STD)
    
    # boundary conditions
    loadingset = domain.getLoadingSet()

    #Physical Line(101) - clamped side of the beam
    loadingset.define(groupset(19), Field1D(TX,RE))
    loadingset.define(groupset(19), Field1D(TY,RE))
    #Physical Line(102) - free surface of the beam  
    #Physical Line(103) - upper surface of the beam (for tests only)


    mim = metafor.getMechanicalIterationManager()
    mim.setMaxNbOfIterations(4)
    mim.setResidualTolerance(1.0e-7)

    ti = AlphaGeneralizedTimeIntegration(metafor)
    metafor.setTimeIntegration(ti)

    # visu
    if 0:
        tsm = metafor.getTimeStepManager()
        tsm.setInitialTime(0.0, 1.0)
        tsm.setNextTime(1.0, 1, 1.0)

    p['exporter'] = None
    return metafor




