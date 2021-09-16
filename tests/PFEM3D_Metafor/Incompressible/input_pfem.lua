Problem = {
    id = "IncompNewtonNoT",
	simulationTime = 0,
	verboseOutput = false,
	
	Mesh = {
		hchar = 0.01,
		alpha = 1.2,
		omega = 0.5,
		gamma = 0.6,
		addOnFS = false,
		deleteFlyingNodes = false,
		laplacianSmoothingBoundaries = false,
		boundingBox = {-0.01,-0.01,0.61,100},
		ignoreGroups = {"SolidBase","Solid"},
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
		ReservoirFixed = true,
		FSInterfaceFixed = false
	},
	
	Solver = {
	    id = "PSPG",
		remesh = false,
		adaptDT = false,
		coeffDTincrease = 1,
		coeffDTDecrease = 1,
		initialDT = 0,
		maxDT = 0,
		
		MomContEq = {
			minRes = 1e-6,
			maxIter = 25,
			bodyForce = {0,-9.81},
			residual = "Ax_f",
			BC = {

			}
		}
	}
}

function Problem.IC:initStates(pos)
	return {0,0,0,0,0}
end

function Problem.Solver.MomContEq.BC:ReservoirV(pos,t)
	return {0,0}
end

function Problem.Solver.MomContEq.BC:FSInterfaceV(pos,t)
	return {0,0}
end

function Problem.Solver.MomContEq.BC:SolidBaseV(pos,t)
	return {0,0}
end