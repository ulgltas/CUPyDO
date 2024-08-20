-- Problem Parameters

Problem = {}
Problem.autoRemeshing = false
Problem.verboseOutput = false
Problem.simulationTime = math.huge
Problem.id = 'WCompNewtonNoT'

-- Mesh Parameters

Problem.Mesh = {}
Problem.Mesh.alpha = 1
Problem.Mesh.omega = 0.7
Problem.Mesh.gamma = 0.5
Problem.Mesh.hchar = 0.01
Problem.Mesh.gammaFS = 0.5
Problem.Mesh.gammaBound = 0.2
Problem.Mesh.minHeightFactor = 1e-2

Problem.Mesh.addOnFS = true
Problem.Mesh.keepFluidElements = true
Problem.Mesh.deleteFlyingNodes = true
Problem.Mesh.deleteBoundElements = {'FSInterface'}
Problem.Mesh.boundingBox = {0, 0, 0.6, 100}
Problem.Mesh.exclusionZones = {}

Problem.Mesh.remeshAlgo = 'CGAL'
Problem.Mesh.mshFile = 'geometryF.msh'

-- Extractor Parameters

Problem.Extractors = {}

Problem.Extractors[0] = {}
Problem.Extractors[0].kind = 'GMSH'
Problem.Extractors[0].writeAs = 'NodesElements'
Problem.Extractors[0].outputFile = 'pfem/fluid.msh'
Problem.Extractors[0].whatToWrite = {'p', 'velocity'}
Problem.Extractors[0].timeBetweenWriting = math.huge

Problem.Extractors[1] = {}
Problem.Extractors[1].kind = 'Global'
Problem.Extractors[1].whatToWrite = 'mass'
Problem.Extractors[1].outputFile = 'mass.txt'
Problem.Extractors[1].timeBetweenWriting = math.huge

-- Material Parameters

Problem.Material = {}
Problem.Material.p0 = 0
Problem.Material.mu = 1e-3
Problem.Material.K0p = 7.6
Problem.Material.K0 = 2.2e+7
Problem.Material.gamma = 7e-2
Problem.Material.rhoStar = 1000

-- Solver Parameters

Problem.Solver = {}
Problem.Solver.id = 'CDS_dpdt'
Problem.Solver.securityCoeff = 0.2

Problem.Solver.adaptDT = true
Problem.Solver.maxDT = math.huge
Problem.Solver.initialDT = math.huge
Problem.Solver.maxRemeshDT = math.huge

-- Momentum Continuity Equation

Problem.Solver.MomEq = {}
Problem.Solver.ContEq = {}
Problem.Solver.MomEq.pExt = 0
Problem.Solver.MomEq.bodyForce = {0, -9.81}
Problem.Solver.ContEq.stabilization = 'CLS'

-- Momentum Continuity BC

Problem.IC = {}
Problem.Solver.MomEq.BC = {}
Problem.Solver.ContEq.BC = {}
Problem.Solver.MomEq.BC['FSInterfaceVExt'] = true

function Problem.IC.initStates(x, y, z)
	return {0, 0, 0, Problem.Material.rhoStar, 0, 0}
end

function Problem.Solver.MomEq.BC.ReservoirV(x, y, z, t)
	return 0, 0
end

function Problem.Solver.MomEq.BC.PolytopeV(x, y, z, t)
	return 0, 0
end