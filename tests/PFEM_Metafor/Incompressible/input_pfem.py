import pfem.tools.link2vtk as v
import pfem
import os

# %% Main Function

def getPfem():
    
    global w
    w = pfem

    file = os.path.dirname(os.path.realpath(__file__))
    file = file+'\geometry.msh'

    pbl = w.Problem()
    msh = w.MshData(pbl)
    msh.load(file)
    
    # Problem Parameters

    toll = 1e-6
    nItMax = 25
    alpha = 1.2
    rho0 = 1000
    mu = 0.001
    g = -9.81
    
    # Assign the Parameters

    pbl.scalingU = 1
    pbl.bodyForceY = g
    pbl.alpha = alpha
    pbl.rho0 = rho0
    pbl.mu = mu

    # Makes Solution Scheme

    solScheme = w.SchemeMonolithicPSPG(msh,pbl)
    
    # Time Integration Scheme
    
    scheme = w.TimeIntegration(msh,pbl,solScheme)
    scheme.addRemoveNodesOption = True
    scheme.savefreq = 10
    scheme.gamma = 0.6
    scheme.omega = 0.5
    
    # Algorithm and Convergence
    
    convCriterion = w.ForceBalanceCriterion(msh,pbl,toll)
    nonLinAlgo = w.PicardAlgorithm(solScheme,convCriterion,nItMax)
    
    # FSI Parameters with 16 = FSInterface

    bndno = 16
    contactTag = w.Tag(100,"Contact",2)
    msh.ntags["Contact"] = contactTag
    msh.ptags[100] = contactTag

    # Medium Parameters
    
    w.Medium(msh,"SSContact",0,0,0,4)
    w.Medium(msh,"FSInterface",0,0,3)
    w.Medium(msh,"Reservoir",mu,rho0,1)
    w.Medium(msh,"Fluid",mu,rho0,1)
    
    # Boundary Conditions

    w.Boundary(msh,"FreeSurface",3,pbl.extP)
    w.Boundary(msh,"Reservoir",1,0)
    w.Boundary(msh,"Reservoir",2,0)
    w.Boundary(msh,"FSInterface",1,0)
    w.Boundary(msh,"FSInterface",2,0)

    # GUI parameters
    
    gui = v.Link2VTK(msh,scheme,False,True)
    
    # Post Processing

    extManager = pfem.ExtractorsManager(msh)
    extManager.add(1,pfem.MassExtractor(msh,pbl,"Fluid"))
    
    # Output of the function
    
    class Module(object):
        def __init__(self):
            
            self.convCriterion = convCriterion
            self.extManager = extManager
            self.contactTag = contactTag
            self.nonLinAlgo = nonLinAlgo
            self.solScheme = solScheme
            self.scheme = scheme
            self.bndno = bndno
            self.msh = msh
            self.pbl = pbl
            self.gui = gui
            self.w = w
        
    return Module()

# %% Dummy Function

def getRealTimeExtractorsList(pfem):
    return []