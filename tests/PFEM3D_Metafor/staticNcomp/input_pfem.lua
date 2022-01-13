Problem = {
    id = "IncompNewtonNoT",
	FSInterface = "FSInterface",
	simulationTime = math.huge,
	verboseOutput = false,
	useCupydo = true,

	Mesh = {
		alpha = 1.2,
		omega = 0.5,
		gamma = 0.6,
		hchar = 0.02,
		addOnFS = false,
		deleteFlyingNodes = false,
		laplacianSmoothingBoundaries = false,
		boundingBox = {-0.01,-1,1.01,1},

		remeshAlgo = "GMSH",
		mshFile = "geometry.msh",
		exclusionGroups = {"FSInterface"},
		ignoreGroups = {"Solid"},
		exclusionZones = {}
	},

	Extractors = {
		{
			kind = "GMSH",
			outputFile = "fluid.msh",
			writeAs = "NodesElements",
			whatToWrite = {"p","velocity"},
			timeBetweenWriting = math.huge
		}
	},

	Material = {
		rho = 1000,
		mu = 1e-3,
		gamma = 0
	},

	IC = {
		WallFixed = true,
		FSInterfaceFixed = false
	},

	Solver = {
	    id = "PSPG",
		coeffDTDecrease = 2,
		coeffDTincrease = 1.5,
		initialDT = math.huge,
		maxDT = math.huge,
		adaptDT = true,
		
		MomContEq = {
			nlAlgo = "Picard",
			residual = "Ax_f",
			bodyForce = {0,-9.81},
			minRes = 1e-8,
			maxIter = 25,
			BC = {}
		}
	}
}

function Problem.IC:initStates(pos)
	return {0,0,0}
end

function Problem.Solver.MomContEq.BC:WallV(pos,t)
	return {0,0}
end

function Problem.Solver.MomContEq.BC:FSInterfaceV(pos,t)
	return nil
end