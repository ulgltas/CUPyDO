# -*- coding: latin-1; -*-
# bends a simple beam from a gmsh file

from wrap import *

metafor = None

def params(q={}):
    """ default model parameters
    """
    p={}
    p['tolNR']      = 1.0e-6        # Newton-Raphson tolerance
    p['tend']       = 1.            # final time
    p['dtmax']      = 0.005          # max time step
    p['bndno']      = 16            # interface boundary number
    
    # BC type
    #p['bctype']     = 'pressure'     # uniform pressure
    #p['bctype']     = 'deadload'     # uniform nodal load
    #p['bctype']     = 'pydeadload1'  # uniform nodal load (python)  
    p['bctype']     = 'pydeadloads'  # variable loads
    #p['bctype']     = 'slave'     # variable loads (mpi)
    
    p['extractor'] = None
                                       
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
    f = os.path.join(os.path.dirname(__file__), "PingPongCylinder_depth_1D.msh")
    importer = GmshImport(f, domain)
    importer.execute2D()

    groupset = domain.getGeometry().getGroupSet()    

    # solid elements / material
    interactionset = domain.getInteractionSet()

    app1 = FieldApplicator(1)
    app1.push( groupset(18) )  # physical group 18: cylinder
    interactionset.add( app1 )
    
    materset = domain.getMaterialSet()
    materset.define( 1, ElastHypoMaterial )
    mater1 = materset(1)
    mater1.put(MASS_DENSITY,    80.0)  # [kg/m�]
    mater1.put(ELASTIC_MODULUS, 100.0e6)  # [Pa]
    mater1.put(POISSON_RATIO,   0.4)   # [-]

    prp = ElementProperties(TriangleVolume2DElement)
    app1.addProperty(prp)
    prp.put (MATERIAL, 1)
    prp.put(GRAVITY_Y, -9.81) # m/s2
    prp.put(CAUCHYMECHVOLINTMETH,VES_CMVIM_STD)
    
    mim = metafor.getMechanicalIterationManager()
    mim.setMaxNbOfIterations(10)
    mim.setResidualTolerance(p['tolNR'])

    ti = AlphaGeneralizedTimeIntegration(metafor)
    metafor.setTimeIntegration(ti)

    # visu
    if 0:
        tsm = metafor.getTimeStepManager()
        tsm.setInitialTime(0.0, 1.0)
        tsm.setNextTime(1.0, 1, 1.0)

    # results
    vmgr = metafor.getValuesManager()
    vmgr.add(1, MiscValueExtractor(metafor, EXT_T), 'time_mtf')
    vmgr.add(2, DbNodalValueExtractor(groupset(18), Field1D(TY,AB), SortByDist0(0.,0.16,0.), 1), 'yPosUp')
    vmgr.add(3, DbNodalValueExtractor(groupset(18), Field1D(TY,AB), SortByDist0(0.,0.14,0.), 1), 'yPosCenter')
    vmgr.add(4, DbNodalValueExtractor(groupset(18), Field1D(TY,RE), SortByDist0(0.,0.16,0.), 1), 'dyPosUp')
    vmgr.add(5, DbNodalValueExtractor(groupset(18), Field1D(TY,RE), SortByDist0(0.,0.14,0.), 1), 'dyPosCenter')
    vmgr.add(6, InteractionGravityCenterAndMassValueExtractor(app1))
    
    return metafor




