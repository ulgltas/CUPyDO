from wrap import FieldApplicator,ElementProperties,Field1D,AlphaGeneralizedTimeIntegration
from wrap import MATERIAL,CAUCHYMECHVOLINTMETH,VES_CMVIM_STD,TX,TY,RE
from wrap import MASS_DENSITY,ELASTIC_MODULUS,POISSON_RATIO
from wrap import ElastHypoMaterial,TriangleVolume2DElement
from toolbox.gmsh import GmshImport
from wrap import Metafor
import os

# %% Physical group 300 = FSInterface

def params(q={}):

    p = {}
    p['tend'] = 0
    p['dtmax'] = 0
    p['bndno'] = 300
    p['tolNR'] = 1e-6
    p['saveAllFacs'] = False
    p['bctype'] = 'pydeadloads'

    p.update(q)
    return p

# %% Main Function

def getMetafor(p={}):
    
    global metafor
    metafor = Metafor()
    p = params(p)
    
    # Problem parameters

    domain = metafor.getDomain()
    geometry = domain.getGeometry()
    geometry.setDimPlaneStrain(1)
    
    # Imports the mesh

    f = os.path.join(os.path.dirname(__file__),"geometry.msh")
    importer = GmshImport(f,domain)
    importer.execute()
    
    # Group and interaction sets

    groupset = domain.getGeometry().getGroupSet()
    interactionset = domain.getInteractionSet()
    
    # Physical group 200 = Solid

    app = FieldApplicator(1)
    app.push(groupset(200))
    interactionset.add(app)
    
    # Material parameters

    materset = domain.getMaterialSet()
    materset.define(1,ElastHypoMaterial)
    mater1 = materset(1)
    mater1.put(MASS_DENSITY,2500)
    mater1.put(ELASTIC_MODULUS,2e4)
    mater1.put(POISSON_RATIO,0)
    
    # FEM properties

    prp = ElementProperties(TriangleVolume2DElement)
    app.addProperty(prp)
    prp.put (MATERIAL,1)
    prp.put(CAUCHYMECHVOLINTMETH,VES_CMVIM_STD)
    
    # Boundary conditions 600 = SolidBase
    
    loadingset = domain.getLoadingSet()
    loadingset.define(groupset(600),Field1D(TX,RE))
    loadingset.define(groupset(600),Field1D(TY,RE))

    # Algorithm parameters

    mim = metafor.getMechanicalIterationManager()
    mim.setMaxNbOfIterations(25)
    mim.setResidualTolerance(p['tolNR'])
    ti = AlphaGeneralizedTimeIntegration(metafor)
    metafor.setTimeIntegration(ti)
    
    return metafor

# %% A Dummy Function

def getRealTimeExtractorsList(mtf):
    
    extractorsList = []
    return extractorsList