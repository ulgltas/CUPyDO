import toolbox.gmsh as gmsh
import wrap as w
import os

# Physical group 1 = FSInterface

def params(parm):

    parm['bndno'] = 1
    return parm

metafor = None
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

    # Defines the solid domain

    app = w.FieldApplicator(1)
    app.push(groups['Solid'])
    interactionset.add(app)
    
    # Material parameters

    materset.define(1,w.ElastHypoMaterial)
    materset(1).put(w.ELASTIC_MODULUS,1e6)
    materset(1).put(w.MASS_DENSITY,2500)
    materset(1).put(w.POISSON_RATIO,0)
    
    # Finite element properties

    prp1 = w.ElementProperties(w.Volume2DElement)
    prp1.put(w.CAUCHYMECHVOLINTMETH,w.VES_CMVIM_STD)
    prp1.put(w.MATERIAL,1)
    app.addProperty(prp1)

    # Elements for surface traction

    prp2 = w.ElementProperties(w.NodStress2DElement)
    load = w.NodInteraction(2)
    load.push(groups['FSInterface'])
    load.addProperty(prp2)
    interactionset.add(load)
    
    # Boundary conditions
    
    loadingset.define(groups['SolidBase'],w.Field1D(w.TX,w.RE))
    loadingset.define(groups['SolidBase'],w.Field1D(w.TY,w.RE))

    # Mechanical time integration

    ti = w.AlphaGeneralizedTimeIntegration(metafor)
    metafor.setTimeIntegration(ti)

    # Mechanical iterations

    mim.setMaxNbOfIterations(4)
    mim.setResidualTolerance(1e-7)

    # Time step iterations

    tscm = w.NbOfMechNRIterationsTimeStepComputationMethod(metafor)
    tsm.setTimeStepComputationMethod(tscm)
    tscm.setTimeStepDivisionFactor(2)
    tscm.setNbOptiIte(25)

    # Parameters for CUPyDO

    parm['interactionM'] = load

    ext = w.GmshExporter(metafor, 'solid')
    ext.add(w.IFNodalValueExtractor(groups['Solid'], w.IF_EVMS))
    ext.add(w.IFNodalValueExtractor(groups['Solid'], w.IF_P))
    parm['exporter'] = ext

    return metafor