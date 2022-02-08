Problem = {
    id = "IncompNewtonNoT",
	FSInterface = "FSInterface",
	simulationTime = math.huge,
	verboseOutput = false,
	autoRemeshing = false,
	
	Mesh = {
		alpha = 1.2,
		omega = 0.5,
		gamma = 0.6,
		hchar = 0.015,
		addOnFS = false,
		deleteFlyingNodes = false,
		laplacianSmoothingBoundaries = false,
		boundingBox = {-0.01,-0.01,0.6,100},

		remeshAlgo = "GMSH",
		mshFile = "geometry.msh",
		exclusionGroups = {"FSInterface"},
		ignoreGroups = {"Solid"},
		exclusionZones = {}
	},
	
	Extractors = {
		{
			kind = "GMSH",
			outputFile = "pfem/fluid.msh",
			writeAs = "NodesElements",
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
		mu = 1e-3,
		gamma = 0
	},

	IC = {
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
			maxIter = 20,
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