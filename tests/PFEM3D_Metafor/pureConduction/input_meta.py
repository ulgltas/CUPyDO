import toolbox.gmsh as gmsh
import numpy as np
import wrap as w
import os

# Physical group 2 = FSInterface

def params(parm):

    parm['bndno'] = 2
    return parm

metafor = None
def getMetafor(parm):

    global metafor
    if metafor: return metafor

    w.StrVectorBase.useTBB()
    w.StrMatrixBase.useTBB()
    w.ContactInteraction.useTBB()
    
    # Group and interaction sets

    metafor = w.Metafor()
    domain = metafor.getDomain()
    tsm = metafor.getTimeStepManager()
    materset = domain.getMaterialSet()
    loadingset = domain.getLoadingSet()
    tim = metafor.getThermalIterationManager()
    solvermanager = metafor.getSolverManager()
    interactionset = domain.getInteractionSet()
    mim = metafor.getMechanicalIterationManager()
    initcondset = metafor.getInitialConditionSet()

    # Dimension and DSS solver

    domain.getGeometry().setDimPlaneStrain(1)
    solvermanager.setSolver(w.DSSolver())
    
    # Imports the mesh

    mshFile = os.path.join(os.path.dirname(__file__), 'geometryS.msh')
    importer = gmsh.GmshImport(mshFile, domain)
    groups = importer.groups
    importer.execute()

    # Defines the solid domain

    app = w.FieldApplicator(1)
    app.push(groups['Solid'])
    interactionset.add(app)

    # Solid material parameters

    materset.define(1, w.TmElastHypoMaterial)
    materset(1).put(w.ELASTIC_MODULUS, 1e7)
    materset(1).put(w.THERM_EXPANSION, 0)
    materset(1).put(w.HEAT_CAPACITY, 1)
    materset(1).put(w.MASS_DENSITY, 8e3)
    materset(1).put(w.POISSON_RATIO, 0)
    materset(1).put(w.CONDUCTIVITY, 10)
    materset(1).put(w.DISSIP_TE, 0)
    materset(1).put(w.DISSIP_TQ, 0)

    # Finite element properties

    prp1 = w.ElementProperties(w.TmVolume2DElement)
    prp1.put(w.CAUCHYMECHVOLINTMETH, w.VES_CMVIM_SRIPR)
    prp1.put(w.STIFFMETHOD, w.STIFF_ANALYTIC)
    prp1.put(w.MATERIAL, 1)
    app.addProperty(prp1)

    # Elements for surface heat flux

    prp2 = w.ElementProperties(w.NodHeatFlux2DElement)
    heat = w.NodInteraction(2)
    heat.push(groups['FSInterface'])
    heat.addProperty(prp2)
    interactionset.add(heat)

    # Top heat flux boundary conditions

    fun = lambda x : np.exp(-np.square(x-0.5)/np.square(0.1)/2)
    F = w.PythonOneParameterFunction(fun)

    prp3 = w.ElementProperties(w.TmHeatFlux2DElement)
    prp3.put(w.HEATEL_VALUE, 1e3)
    prp3.depend(w.HEATEL_VALUE, F, w.Field1D(w.TX, w.AB))

    inter = w.HeatInteraction(3)
    inter.push(groups['Top'])
    inter.addProperty(prp3)
    interactionset.add(inter)

    # Other boundary conditions

    initcondset.define(groups['Solid'], w.Field1D(w.TO, w.AB), 300)
    loadingset.define(groups['FSInterface'], w.Field1D(w.TX, w.RE))
    loadingset.define(groups['FSInterface'], w.Field1D(w.TY, w.RE))

    # Mechanical and thermal time integration

    ti_M = w.AlphaGeneralizedTimeIntegration(metafor)
    ti_T = w.TrapezoidalThermalTimeIntegration(metafor)

    ti = w.StaggeredTmTimeIntegration(metafor)
    ti.setMechanicalTimeIntegration(ti_M)
    ti.setThermalTimeIntegration(ti_T) 
    metafor.setTimeIntegration(ti)

    # Mechanical and thermal iterations

    mim.setMaxNbOfIterations(25)
    mim.setResidualTolerance(1e-6)

    tim.setMaxNbOfIterations(25)
    tim.setResidualTolerance(1e-6)

    # Time step iterations
    
    tscm = w.NbOfStaggeredTmNRIterationsTimeStepComputationMethod(metafor)
    tsm.setTimeStepComputationMethod(tscm)
    tscm.setTimeStepDivisionFactor(2)
    tscm.setNbOptiIte(25)

    # Parameters for CUPyDO

    parm['interactionT'] = heat
    parm['FSInterface'] = groups['FSInterface']

    ext = w.GmshExporter(metafor, 'solid')
    ext.add(w.DbNodalValueExtractor(groups['Solid'], w.Field1D(w.TO, w.AB)))
    ext.add(w.DbNodalValueExtractor(groups['Solid'], w.Field1D(w.TO, w.RE)))

    parm['exporter'] = ext

    return metafor