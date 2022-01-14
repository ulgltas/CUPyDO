Problem = {
    id = "WCompNewtonNoT",
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
		boundingBox = {-0.01,-0.01,0.61,100},

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
		rhoStar = 1000,
        K0 = 2.2e+7,
        K0p = 7.6,
		mu = 1e-3,
		gamma = 0
	},
	
	IC = {
		ReservoirFixed = true,
		FSInterfaceFixed = false,
	},
	
	Solver = {
	    id = "CDS_dpdt",
		securityCoeff = 0.2,
		initialDT = math.huge,
		maxDT = math.huge,
		adaptDT = true,

        MomEq = {
            bodyForce = {0,-9.81},
            BC = {}
        },
        
        ContEq = {
            stabilization = "Meduri",
            BC = {}
		}
	}
}

function Problem.IC:initStates(pos)
	return {0,0,0,Problem.Material.rhoStar,0,0}
end

function Problem.Solver.MomEq.BC:ReservoirV(pos,t)
	return {0,0}
end

function Problem.Solver.MomEq.BC:FSInterfaceV(pos,t)
	return nil
end