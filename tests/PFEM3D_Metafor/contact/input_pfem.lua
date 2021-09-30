Problem = {
    id = "IncompNewtonNoT",
	simulationTime = 0,
	verboseOutput = false,
	
	Mesh = {
		hchar = 0.005,
		alpha = 1.2,
		omega = 0.5,
		gamma = 0.6,
		addOnFS = false,
		deleteFlyingNodes = false,
		laplacianSmoothingBoundaries = false,
		boundingBox = {-0.111,-0.151,0.111,0.211},
		ignoreGroups = {"SolidBaseR","SolidBaseL","SolidR","SolidL","Ball"},
		exclusionZones = {},
		mshFile = "geometry.msh"
	},
	
	Extractors = {
		{
			kind = "GMSH",
			outputFile = "results.msh",
			timeBetweenWriting = 0.01,
			whatToWrite = {"p","velocity"},
			writeAs = "NodesElements" 
		}
	},
	
	Material = {
		mu = 1e-3,
		rho = 1000,
		gamma = 0
	},
	
	IC = {
		InletFixed = true,
		OutletFixed = true,
		ReservoirFixed = false,
		FSInterfaceFixed = false,

		FSInterfaceRFixed = false,
		FSInterfaceBFixed = false
	},

	Solver = {
	    id = "PSPG",
		adaptDT = false,
		coeffDTincrease = 1,
		coeffDTDecrease = 1,
		initialDT = 0,
		maxDT = 0,
		
		MomContEq = {
			minRes = 1e-6,
			maxIter = 25,
			bodyForce = {0,0},
			residual = "Ax_f",
			BC = {

			}
		}
	}
}

function Problem.IC:initStates(pos)
	return {0,0,0}
end

function Problem.Solver.MomContEq.BC:ReservoirV(pos,t)
	return {0,0}
end

function Problem.Solver.MomContEq.BC:FSInterfaceV(pos,t)
	return {0,0}
end

function Problem.Solver.MomContEq.BC:FSInterfaceRV(pos,t)
	return {0,0}
end

function Problem.Solver.MomContEq.BC:FSInterfaceBV(pos,t)
	return {0,0}
end

function Problem.Solver.MomContEq.BC:InletV(pos,t)
	return {0,-0.01}
end

function Problem.Solver.MomContEq.BC:OutletV(pos,t)
	return {0,-0.01}
end