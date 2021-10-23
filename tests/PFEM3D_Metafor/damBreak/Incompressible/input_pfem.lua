Problem = {
    id = "IncompNewtonNoT",
	simulationTime = math.huge,
	verboseOutput = false,
	useCupydo = true,
	
	Mesh = {
		alpha = 1.2,
		omega = 0.5,
		gamma = 0.6,
		hchar = 0.01,
		addOnFS = false,
		deleteFlyingNodes = false,
		laplacianSmoothingBoundaries = false,
		boundingBox = {-0.01,-0.01,0.61,100},
		ignoreGroups = {"SolidBase","Solid"},
		mshFile = "geometry.msh",
		exclusionZones = {}
	},
	
	Extractors = {
		{
			kind = "GMSH",
			writeAs = "NodesElements",
			timeBetweenWriting = 0.01,
			whatToWrite = {"p","velocity"},
			outputFile = "fluid.msh"
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
		adaptDT = true,
		coeffDTincrease = 1.5,
		coeffDTDecrease = 2,
		initialDT = 0.1,
		maxDT = 0.1,
		
		MomContEq = {
			maxIter = 25,
			minRes = 1e-8,
			bodyForce = {0,-9.81},
			residual = "Ax_f",
			nlAlgo = "Picard",
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