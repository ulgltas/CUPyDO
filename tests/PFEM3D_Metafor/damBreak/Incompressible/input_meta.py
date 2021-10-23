from toolbox.gmsh import GmshImport
import wrap as w
import os

# %% Physical group 16 = FSInterface

def params(input):

    input['bndno'] = 16
    input['saveAllFacs'] = False
    input['bctype'] = 'pydeadloads'
    return input

# %% Parallel Computing

metafor = None
w.StrVectorBase.useTBB()
w.StrMatrixBase.useTBB()
w.ContactInteraction.useTBB()

# %% Main Function

def getMetafor(input):

    global metafor
    if metafor: return metafor
    
    # Group and interaction sets

    metafor = w.Metafor()
    domain = metafor.getDomain()
    materset = domain.getMaterialSet()
    loadingset = domain.getLoadingSet()
    solvermanager = metafor.getSolverManager()
    interactionset = domain.getInteractionSet()
    mim = metafor.getMechanicalIterationManager()

    # Plane strain and DSS solver

    domain.getGeometry().setDimPlaneStrain(1)
    solvermanager.setSolver(w.DSSolver())
    
    # Imports the mesh

    mshFile = os.path.join(os.path.dirname(__file__),"geometry.msh")
    importer = GmshImport(mshFile,domain)
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

    prp = w.ElementProperties(w.Volume2DElement)
    prp.put(w.CAUCHYMECHVOLINTMETH,w.VES_CMVIM_STD)
    prp.put(w.MATERIAL,1)
    app.addProperty(prp)
    
    # Boundary conditions
    
    loadingset.define(groups['SolidBase'],w.Field1D(w.TX,w.RE))
    loadingset.define(groups['SolidBase'],w.Field1D(w.TY,w.RE))

    # Algorithm parameters

    mim.setMaxNbOfIterations(25)
    mim.setResidualTolerance(1e-8)
    ti = w.AlphaGeneralizedTimeIntegration(metafor)
    metafor.setTimeIntegration(ti)

    return metafor

# %% A Dummy Function

def getRealTimeExtractorsList(mtf):
    
    extractorsList = []
    return extractorsList