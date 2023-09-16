-- Problem Parameters

Problem = {}
Problem.verboseOutput = true
Problem.autoRemeshing = false
Problem.simulationTime = math.huge
Problem.id = 'Conduction'

-- Mesh Parameters

Problem.Mesh = {}
Problem.Mesh.alpha = 1.2
Problem.Mesh.omega = 0.5
Problem.Mesh.gamma = 0.6
Problem.Mesh.hchar = 0.01
Problem.Mesh.gammaFS = 0.2
Problem.Mesh.addOnFS = false
Problem.Mesh.minAspectRatio = 1e-3
Problem.Mesh.keepFluidElements = true
Problem.Mesh.deleteFlyingNodes = false
Problem.Mesh.deleteBoundElements = false
Problem.Mesh.laplacianSmoothingBoundaries = false
Problem.Mesh.boundingBox = {0,-10,1,1}
Problem.Mesh.exclusionZones = {}

Problem.Mesh.remeshAlgo = 'GMSH'
Problem.Mesh.mshFile = 'geometryF.msh'
Problem.Mesh.exclusionGroups = {}
Problem.Mesh.ignoreGroups = {}

-- Extractor Parameters

Problem.Extractors = {}

Problem.Extractors[0] = {}
Problem.Extractors[0].kind = 'GMSH'
Problem.Extractors[0].writeAs = 'NodesElements'
Problem.Extractors[0].outputFile = 'pfem/output.msh'
Problem.Extractors[0].whatToWrite = {'T'}
Problem.Extractors[0].timeBetweenWriting = math.huge

Problem.Extractors[1] = {}
Problem.Extractors[1].kind = 'Global'
Problem.Extractors[1].whatToWrite = 'mass'
Problem.Extractors[1].outputFile = 'mass.txt'
Problem.Extractors[1].timeBetweenWriting = math.huge

-- Material Parameters

Problem.Material = {}
Problem.Material.mu = 1e-3
Problem.Material.gamma = 0
Problem.Material.rho = 1000
Problem.Material.epsRad = 0
Problem.Material.sigmaRad = 5.670374419e-8
Problem.Material.R = 8.31446261815324
Problem.Material.alphaLin = 69e-6
Problem.Material.DgammaDT = 0
Problem.Material.Tinf = 300
Problem.Material.DmuDT = 0
Problem.Material.DcpDT = 0
Problem.Material.DkDT = 0
Problem.Material.Tr = 650
Problem.Material.k = 0.6
Problem.Material.cp = 1
Problem.Material.h = 5

-- Solver Parameters

Problem.Solver = {}
Problem.Solver.id = 'PSPG'

Problem.Solver.adaptDT = true
Problem.Solver.maxDT = math.huge
Problem.Solver.initialDT = math.huge
Problem.Solver.coeffDTDecrease = math.huge
Problem.Solver.coeffDTincrease = math.huge
Problem.Solver.solveHeatFirst = true

-- Heat Equation

Problem.Solver.HeatEq = {}
Problem.Solver.HeatEq.residual = 'Ax_f'
Problem.Solver.HeatEq.nlAlgo = 'Picard'
Problem.Solver.HeatEq.sparseSolver = 'CG'

Problem.Solver.HeatEq.maxIter = 25
Problem.Solver.HeatEq.minRes = 1e-8
Problem.Solver.HeatEq.cgTolerance = 1e-16

-- Heat Momentum Continuity BC

Problem.IC = {}
Problem.Solver.HeatEq.BC = {}
Problem.Solver.HeatEq.BC['FSInterfaceTExt'] = true

function Problem.IC.initStates(x,y,z)
	return {300}
end

function Problem.Solver.HeatEq.BC.WallQ(x,y,z,t) 
    return 0,0
end