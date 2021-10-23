Problem = {
    id = "IncompNewtonNoT",
	simulationTime = math.huge,
	verboseOutput = false,
	useCupydo = true,
	
	Mesh = {
		alpha = 1.2,
		omega = 0.5,
		gamma = 0.6,
		hchar = 6e-3,
		addOnFS = false,
		deleteFlyingNodes = false,
		laplacianSmoothingBoundaries = false,
		boundingBox = {-0.101,-0.001,0.101,0.151},
		ignoreGroups = {"SolidBase","Solid"},
		mshFile = "geometry.msh",
		exclusionZones = {}
	},
	
	Extractors = {
		{
			kind = "GMSH",
			outputFile = "fluid.msh",
			timeBetweenWriting = 0.1,
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
		ReservoirFixed = false,
		FSInterfaceFixed = false,
		InletFixed = true
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
			bodyForce = {0,0},
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

function Problem.Solver.MomContEq.BC:InletV(pos,t)
	return {0.1,0}
end