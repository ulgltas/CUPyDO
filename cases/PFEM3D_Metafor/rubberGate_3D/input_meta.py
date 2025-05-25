import toolbox.gmsh as gmsh
import wrap as w
import os

metafor = None

def getMetafor(parm):

    global metafor
    if metafor: return metafor
    metafor = w.Metafor()

    w.StrVectorBase.useTBB()
    w.StrMatrixBase.useTBB()
    w.ContactInteraction.useTBB()

    domain = metafor.getDomain()
    domain.getGeometry().setDim3D()
    metafor.getSolverManager().setSolver(w.DSSolver())

    mshFile = os.path.join(os.path.dirname(__file__), 'geometryS.msh')
    importer = gmsh.GmshImport(mshFile, domain)
    groups = importer.groups
    importer.execute()

    parm['FSI'] = groups['FSI']

    iset = domain.getInteractionSet()
    app = w.FieldApplicator(1)
    app.push(groups['Solid'])
    iset.add(app)

    K = 1.3e6
    G = 2.4e6
    C1 = -1.2e6
    C2 = G/2.0-C1

    materset = domain.getMaterialSet()
    materset.define(1, w.MooneyRivlinHyperMaterial)
    materset(1).put(w.MASS_DENSITY, 1100)
    materset(1).put(w.RUBBER_PENAL, K)
    materset(1).put(w.RUBBER_C1, C1)
    materset(1).put(w.RUBBER_C2, C2)
    
    prp1 = w.ElementProperties(w.Volume3DElement)
    prp1.put(w.CAUCHYMECHVOLINTMETH, w.VES_CMVIM_SRIPR)
    prp1.put(w.STIFFMETHOD, w.STIFF_ANALYTIC)
    prp1.put(w.GRAVITY_Z, -9.81)
    prp1.put(w.MATERIAL, 1)
    app.addProperty(prp1)

    prp2 = w.ElementProperties(w.NodStress3DElement)
    load = w.NodInteraction(2)
    load.push(groups['FSI'])
    load.addProperty(prp2)
    iset.add(load)

    parm['interactionM'] = load
    
    loadset = domain.getLoadingSet()
    loadset.define(groups['Base'], w.Field1D(w.TX, w.RE))
    loadset.define(groups['Base'], w.Field1D(w.TY, w.RE))
    loadset.define(groups['Base'], w.Field1D(w.TZ, w.RE))

    ti = w.AlphaGeneralizedTimeIntegration(metafor)
    metafor.setTimeIntegration(ti)

    mim = metafor.getMechanicalIterationManager()
    mim.setResidualTolerance(1e-6)
    mim.setMaxNbOfIterations(25)

    tsm = metafor.getTimeStepManager()
    tscm = w.NbOfMechNRIterationsTimeStepComputationMethod(metafor)
    tsm.setTimeStepComputationMethod(tscm)
    tscm.setTimeStepDivisionFactor(2)
    tscm.setNbOptiIte(25)

    ext = w.GmshExporter(metafor, 'metafor/output')
    ext.add(w.IFNodalValueExtractor(groups['Solid'], w.IF_P))
    ext.add(w.IFNodalValueExtractor(groups['Solid'], w.IF_EVMS))
    parm['exporter'] = ext

    # Build domain and folder

    domain.build()
    os.makedirs('metafor')
    return metafor