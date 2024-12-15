import toolbox.gmsh as gmsh
import wrap as w
import os

# Physical group 2 = FSI

def params(parm):

    parm['bndno'] = 2
    return parm

metafor = None
def getMetafor(parm):
    '''
    Initialize and return the metafor object
    '''

    global metafor
    if metafor: return metafor
    metafor = w.Metafor()
    
    # Group and interaction sets

    domain = metafor.getDomain()
    iset = domain.getInteractionSet()
    tsm = metafor.getTimeStepManager()
    materset = domain.getMaterialSet()
    loadingset = domain.getLoadingSet()
    solvermanager = metafor.getSolverManager()
    interactionset = domain.getInteractionSet()
    mim = metafor.getMechanicalIterationManager()

    # Dimension and DSS solver

    domain.getGeometry().setDim3D()
    solvermanager.setSolver(w.DSSolver())

    # Imports the mesh

    mshFile = os.path.join(os.path.dirname(__file__), 'geometryS.msh')
    importer = gmsh.GmshImport(mshFile, domain)
    groups = importer.groups
    importer.execute()

    # Defines the solid domain

    app = w.FieldApplicator(1)
    app.push(groups['Solid'])
    iset.add(app)
    
    # Material parameters

    K = 1.3e6
    G = 2.4e6
    C1 = -1.2e6
    C2 = G/2.0-C1

    materset.define(1, w.MooneyRivlinHyperMaterial)
    materset(1).put(w.MASS_DENSITY, 1100)
    materset(1).put(w.RUBBER_PENAL,  K)
    materset(1).put(w.RUBBER_C1, C1)
    materset(1).put(w.RUBBER_C2, C2)

    # Finite element properties

    prp1 = w.ElementProperties(w.Volume3DElement)
    prp1.put(w.CAUCHYMECHVOLINTMETH, w.VES_CMVIM_STD)
    prp1.put(w.STIFFMETHOD, w.STIFF_NUMERIC)
    prp1.put(w.GRAVITY_Z, -9.81)
    prp1.put(w.MATERIAL, 1)
    app.addProperty(prp1)

    # Elements for surface traction

    prp2 = w.ElementProperties(w.NodStress3DElement)
    load = w.NodInteraction(2)
    load.push(groups['FSI'])
    load.addProperty(prp2)
    iset.add(load)

    parm['interactionM'] = load
    
    # Boundary conditions
    
    loadingset.define(groups['Base'], w.Field1D(w.TX, w.RE))
    loadingset.define(groups['Base'], w.Field1D(w.TY, w.RE))
    loadingset.define(groups['Base'], w.Field1D(w.TZ, w.RE))
    loadingset.define(groups['Side'], w.Field1D(w.TY, w.RE))

    # Mechanical time integration

    ti = w.AlphaGeneralizedTimeIntegration(metafor)
    metafor.setTimeIntegration(ti)

    # Mechanical iterations

    mim.setMaxNbOfIterations(25)
    mim.setResidualTolerance(1e-8)

    # Time step iterations

    tscm = w.NbOfMechNRIterationsTimeStepComputationMethod(metafor)
    tsm.setTimeStepComputationMethod(tscm)
    tscm.setTimeStepDivisionFactor(2)
    tscm.setNbOptiIte(25)

    # Nodal Gmsh exporter

    ext = w.GmshExporter(metafor, 'metafor/output')
    ext.add(w.IFNodalValueExtractor(groups['Solid'], w.IF_P))
    ext.add(w.IFNodalValueExtractor(groups['Solid'], w.IF_EVMS))
    parm['exporter'] = ext

    # Build domain and folder

    domain.build()
    os.makedirs('metafor')
    return metafor