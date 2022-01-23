Problem = {
    id = "IncompNewtonNoT",
	FSInterface = "FSInterface",
	simulationTime = math.huge,
	verboseOutput = false,
	useCupydo = true,
	
	Mesh = {
		alpha = 1.2,
		omega = 0.4,
		gamma = 0.8,
		hchar = 0.1,
		addOnFS = false,
		deleteFlyingNodes = false,
		laplacianSmoothingBoundaries = false,
		boundingBox = {-4,-4,4,6.3},

		remeshAlgo = "GMSH",
		mshFile = "geometry.msh",
		exclusionGroups = {"FSInterface","PolyL","PolyR"},
		ignoreGroups = {"Solid"},
		exclusionZones = {}
	},
	
	Extractors = {
		{
			kind = "GMSH",
			writeAs = "NodesElements",
			outputFile = "pfem/fluid.msh",
			whatToWrite = {"p","velocity"},
			timeBetweenWriting = math.huge
		},
        {
            kind = "Global",
            whatToWrite = "mass",
            outputFile = "mass.txt",
            timeBetweenWriting = math.huge
        }
	},

	Material = {
		rho = 1000,
		gamma = 0,
		mu = 100
	},

	IC = {
		PolyLFixed = true,
		PolyRFixed = true,
		ReservoirFixed = true,
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
			minRes = 1e-6,
			maxIter = 25,
			BC = {}
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
	return nil
end