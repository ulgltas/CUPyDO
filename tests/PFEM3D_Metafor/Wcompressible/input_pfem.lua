Problem = {
    id = "WCompNewtonNoT",
	simulationTime = 0,
	verboseOutput = false,
	
	Mesh = {
        hchar = 0.0048,
        alpha = 1.2,
        omega = 0.7,
        gamma = 0.7,
		addOnFS = true,
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
		rhoStar = 1000,
        K0 = 2.2e+7,
        K0p = 7.6,
		gamma = 0
	},
	
	IC = {
		ReservoirFixed = true,
		FSInterfaceFixed = false
	},
	
	Solver = {
	    id = "CDS_dpdt",
		remesh = false,
		adaptDT = false,
		securityCoeff = 0,
		initialDT = 0,
		maxDT = 0,

        MomEq = {
            bodyForce = {0,-9.81},
            BC = {

            }
        },
        
        ContEq = {
            stabilization = "Meduri",
            BC = {

            }
		}
	}
}

function Problem.IC:initStates(pos)

	local rho = Problem.Material.rhoStar
	return {0,0,0,rho,0,0}
end

function Problem.Solver.MomEq.BC:ReservoirV(pos,t)
	return {0,0}
end

function Problem.Solver.MomEq.BC:FSInterfaceV(pos,t)
	return {0,0}
end

function Problem.Solver.MomEq.BC:SolidBaseV(pos,t)
	return {0,0}
end