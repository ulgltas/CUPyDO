from toolbox.gmsh import GmshImport
import wrap as w
import os

# %% Physical group 300 = FSInterface

def params(param={}):

    param['tend'] = 0
    param['bndno'] = 300
    param['dtmax'] = 0.01
    param['tolNR'] = 1e-6
    param['saveAllFacs'] = False
    param['bctype'] = 'pydeadloads'
    return param

# %% Main Function

def getMetafor(param={}):
    
    global metafor
    metafor = w.Metafor()
    param = params(param)
    
    # Group and interaction sets

    domain = metafor.getDomain()
    geometry = domain.getGeometry()
    groupset = geometry .getGroupSet()
    materset = domain.getMaterialSet()
    loadingset = domain.getLoadingSet()
    interactionset = domain.getInteractionSet()
    mim = metafor.getMechanicalIterationManager()
    geometry.setDimPlaneStrain(1)
    
    # Imports the mesh

    f = os.path.join(os.path.dirname(__file__),"geometry.msh")
    importer = GmshImport(f,domain)
    importer.execute()
    
    # Physical group 200 = Solid

    app = w.FieldApplicator(1)
    app.push(groupset(200))
    interactionset.add(app)
    
    # Material parameters

    materset.define(1,w.ElastHypoMaterial)
    mater1 = materset(1)
    mater1.put(w.MASS_DENSITY,2500)
    mater1.put(w.ELASTIC_MODULUS,1e4)
    mater1.put(w.POISSON_RATIO,0)
    
    # FEM properties

    prp = w.ElementProperties(w.TriangleVolume2DElement)
    app.addProperty(prp)
    prp.put(w.MATERIAL,1)
    prp.put(w.CAUCHYMECHVOLINTMETH,w.VES_CMVIM_STD)
    
    # Boundary conditions 600 = SolidBase
    
    loadingset.define(groupset(600),w.Field1D(w.TX,w.RE))
    loadingset.define(groupset(600),w.Field1D(w.TY,w.RE))

    # Algorithm parameters

    mim.setMaxNbOfIterations(25)
    mim.setResidualTolerance(param['tolNR'])
    ti = w.AlphaGeneralizedTimeIntegration(metafor)
    metafor.setTimeIntegration(ti)
    
    return metafor

# %% A Dummy Function

def getRealTimeExtractorsList(mtf):
    
    extractorsList = []
    return extractorsList