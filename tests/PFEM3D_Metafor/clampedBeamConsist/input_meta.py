import toolbox.gmsh as gmsh
import wrap as w
import os

# Physical group 2 = FSInterface

def params(parm):

    parm['bndno'] = 2
    return parm

# Parallel Computing

metafor = None
w.StrVectorBase.useTBB()
w.StrMatrixBase.useTBB()
w.ContactInteraction.useTBB()

# Main Function

def getMetafor(parm):

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

    # Plane strain and DSS solver

    domain.getGeometry().setDimPlaneStrain(1)
    solvermanager.setSolver(w.DSSolver())
    
    # Imports the mesh

    mshFile = os.path.join(os.path.dirname(__file__),"geometryS.msh")
    importer = gmsh.GmshImport(mshFile,domain)
    groups = importer.groups
    importer.execute()

    # Defines the ball domain

    app = w.FieldApplicator(1)
    app.push(groups['Solid'])
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

    # Elements for surface traction

    prp2 = w.ElementProperties(w.NodStress2DElement)
    load = w.NodInteraction(2)
    load.push(groups['FSInterface'])
    load.addProperty(prp2)
    interactionset.add(load)

    # Boundary conditions
    
    loadingset.define(groups['Clamped'],w.Field1D(w.TX,w.RE))
    loadingset.define(groups['Clamped'],w.Field1D(w.TY,w.RE))

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

    # Parameters for CUPyDO

    parm['interacM'] = load
    parm['exporter'] = gmsh.GmshExport('solid.msh',metafor)
    parm['exporter'].addInternalField([w.IF_EVMS,w.IF_P])
    return metafor