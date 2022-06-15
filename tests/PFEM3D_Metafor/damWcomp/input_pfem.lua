-- Problem Parameters

Problem = {}
Problem.autoRemeshing = false
Problem.verboseOutput = false
Problem.simulationTime = math.huge
Problem.id = 'WCompNewtonNoT'

-- FSPC Parameters

Problem.interface = 'FSInterface'
Problem.maxFactor = 1000

-- Mesh Parameters

Problem.Mesh = {}
Problem.Mesh.alpha = 1.2
Problem.Mesh.omega = 0.5
Problem.Mesh.gamma = 0.6
Problem.Mesh.hchar = 0.015
Problem.Mesh.addOnFS = false
Problem.Mesh.keepFluidElements = true
Problem.Mesh.deleteFlyingNodes = false
Problem.Mesh.deleteBoundElements = false
Problem.Mesh.laplacianSmoothingBoundaries = false
Problem.Mesh.boundingBox = {0,0,0.6,0.6}
Problem.Mesh.exclusionZones = {}

Problem.Mesh.remeshAlgo = 'CGAL'
Problem.Mesh.mshFile = 'geometry.msh'
Problem.Mesh.exclusionGroups = {'FSInterface'}
Problem.Mesh.ignoreGroups = {'Solid'}

-- Extractor Parameters

Problem.Extractors = {}

Problem.Extractors[0] = {}
Problem.Extractors[0].kind = 'GMSH'
Problem.Extractors[0].writeAs = 'NodesElements'
Problem.Extractors[0].outputFile = 'fluid.msh'
Problem.Extractors[0].whatToWrite = {'p','velocity'}
Problem.Extractors[0].timeBetweenWriting = math.huge

Problem.Extractors[1] = {}
Problem.Extractors[1].kind = 'Global'
Problem.Extractors[1].whatToWrite = 'mass'
Problem.Extractors[1].outputFile = 'mass.txt'
Problem.Extractors[1].timeBetweenWriting = math.huge

-- Material Parameters

Problem.Material = {}
Problem.Material.mu = 1e-3
Problem.Material.K0p = 7.6
Problem.Material.gamma = 0
Problem.Material.K0 = 2.2e+7
Problem.Material.rhoStar = 1000

-- Initial Conditions

Problem.IC = {}
Problem.IC.ReservoirFixed = true
Problem.IC.FSInterfaceFixed = false

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
Problem.Solver.MomEq.bodyForce = {0,-9.81}
Problem.Solver.ContEq.stabilization = 'Meduri'

-- Momentum Continuity BC

Problem.Solver.MomEq.BC = {}
Problem.Solver.ContEq.BC = {}
Problem.Solver.MomEq.BC['FSInterfaceVExt'] = true

function Problem.IC:initStates(pos)
	return {0,0,0,Problem.Material.rhoStar,0,0}
end

function Problem.Solver.MomEq.BC:ReservoirV(pos,t)
	return {0,0}
end