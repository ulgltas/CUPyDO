import toolbox.gmsh as gmsh
import wrap as w
import os

metafor = None
def getMetafor(parm):
    '''
    Initialize and return the metafor object
    '''

    global metafor
    if metafor: return metafor
    metafor = w.Metafor()

    w.StrVectorBase.useTBB()
    w.StrMatrixBase.useTBB()
    w.ContactInteraction.useTBB()

    # Dimension and DSS solver

    domain = metafor.getDomain()
    domain.getGeometry().setDimPlaneStrain(1)
    metafor.getSolverManager().setSolver(w.DSSolver())
    
    # Imports the mesh

    mshFile = os.path.join(os.path.dirname(__file__), 'geometryS.msh')
    importer = gmsh.GmshImport(mshFile, domain)
    groups = importer.groups
    importer.execute()

    parm['FSI'] = groups['FSI']

    # Defines the solid domain

    iset = domain.getInteractionSet()
    app = w.FieldApplicator(1)
    app.push(groups['Solid'])
    iset.add(app)

    # Material parameters

    matlawset = domain.getMaterialLawSet()
    matlawset.define(1, w.LinearIsotropicHardening)
    matlawset(1).put(w.IH_SIGEL, 5e4)
    matlawset(1).put(w.IH_H, 2e4)

    materset = domain.getMaterialSet()
    materset.define(1, w.EvpIsoHHypoMaterial)
    materset(1).put(w.ELASTIC_MODULUS, 1e6)
    materset(1).put(w.MASS_DENSITY, 2500)
    materset(1).put(w.POISSON_RATIO, 0)
    materset(1).put(w.YIELD_NUM, 1)
    
    # Finite element properties

    prp1 = w.ElementProperties(w.Volume2DElement)
    prp1.put(w.CAUCHYMECHVOLINTMETH, w.VES_CMVIM_SRIPR)
    prp1.put(w.STIFFMETHOD, w.STIFF_ANALYTIC)
    prp1.put(w.MATERIAL, 1)
    app.addProperty(prp1)

    # Elements for surface traction

    prp2 = w.ElementProperties(w.NodStress2DElement)
    load = w.NodInteraction(2)
    load.push(groups['FSI'])
    load.push(groups['Base'])
    load.addProperty(prp2)
    iset.add(load)

    parm['interactionM'] = load
    
    # Boundary conditions
    
    loadset = domain.getLoadingSet()
    loadset.define(groups['Base'], w.Field1D(w.TX, w.RE))
    loadset.define(groups['Base'], w.Field1D(w.TY, w.RE))

    # Mechanical time integration

    ti = w.AlphaGeneralizedTimeIntegration(metafor)
    metafor.setTimeIntegration(ti)

    # Mechanical iterations

    mim = metafor.getMechanicalIterationManager()
    mim.setResidualTolerance(1e-7)
    mim.setMaxNbOfIterations(25)

    # Time step iterations

    tsm = metafor.getTimeStepManager()
    tscm = w.NbOfMechNRIterationsTimeStepComputationMethod(metafor)
    tsm.setTimeStepComputationMethod(tscm)
    tscm.setTimeStepDivisionFactor(2)
    tscm.setNbOptiIte(25)

    # Nodal Gmsh exporter

    ext = w.GmshExporter(metafor, 'metafor/output')
    ext.add(w.IFNodalValueExtractor(groups['Solid'], w.IF_EVMS))
    ext.add(w.DbNodalValueExtractor(groups['Solid'], w.Field1D(w.TX, w.GF1)))
    parm['exporter'] = ext

    # Build domain and folder

    domain.build()
    os.makedirs('metafor')
    return metafor