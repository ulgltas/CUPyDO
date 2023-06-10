import toolbox.gmshOld as gmshOld
import toolbox.gmsh as gmsh
import wrap as w
import os

# %% Physical group 3 = FSInterface

def params(p):

    p['bndno'] = 3
    p['saveAllFacs'] = False
    p['bctype'] = 'pydeadloads'
    p['exporter'] = gmsh.GmshExport('metafor/solid.msh',metafor)
    p['exporter'].addInternalField([w.IF_EVMS,w.IF_P])
    return p

# %% Parallel Computing

metafor = None
w.StrVectorBase.useTBB()
w.StrMatrixBase.useTBB()
w.ContactInteraction.useTBB()

# %% Main Function

def getMetafor(p):

    global metafor
    if metafor: return metafor
    
    # Group and interaction sets

    metafor = w.Metafor()
    domain = metafor.getDomain()
    tsm = metafor.getTimeStepManager()
    materset = domain.getMaterialSet()
    loadingset = domain.getLoadingSet()
    solvermanager = metafor.getSolverManager()
    interactionset = domain.getInteractionSet()
    mim = metafor.getMechanicalIterationManager()
    groupset = domain.getGeometry().getGroupSet()

    # Plane strain and DSS solver

    domain.getGeometry().setDimPlaneStrain(1)
    solvermanager.setSolver(w.DSSolver())
    
    # Imports the mesh

    mshFile = os.path.join(os.path.dirname(__file__),"beam.msh")
    importer = gmshOld.GmshImport(mshFile,domain)
    importer.execute2D()

    # Defines the ball domain

    app = w.FieldApplicator(1)
    app.push(groupset(1))
    interactionset.add(app)

    # Solid material parameters

    materset.define(1,w.ElastHypoMaterial)
    materset(1).put(w.ELASTIC_MODULUS,1e7)
    materset(1).put(w.MASS_DENSITY,8e3)
    materset(1).put(w.POISSON_RATIO,0)

    # Finite element properties

    prp = w.ElementProperties(w.Volume2DElement)
    prp.put(w.CAUCHYMECHVOLINTMETH,w.VES_CMVIM_STD)
    prp.put(w.MATERIAL,1)
    app.addProperty(prp)

    # Boundary conditions
    
    loadingset.define(groupset(5),w.Field1D(w.TX,w.RE))
    loadingset.define(groupset(5),w.Field1D(w.TY,w.RE))

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

    return metafor