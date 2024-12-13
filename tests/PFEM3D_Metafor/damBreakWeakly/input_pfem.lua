Problem = {

    -- Main simulation parameters

    mechanical = true,
    verboseOutput = false,
    autoRemeshing = false,
    simulationTime = math.huge
}

Problem.Mesh = {

    -- Input mesh and bounding box

    remeshAlgo = 'CGAL',
    mshFile = 'geometryF.msh',
    deleteBoundElements = {'FSI'},
    boundingBox = {0, 0, 0.6, 100},
    exclusionZones = {},

    -- Remeshing internal parameters

    alpha = 1.2,
    omega = 0.5,
    gamma = 0.4,
    hchar = 0.01,
    gammaFS = 0.4,
    gammaEdge = 0.2,
    minHeightFactor = 1e-3,

    -- Enable or disable algorithms

    addOnFS = true,
    keepFluidElements = true,
    deleteFlyingNodes = true
}

Problem.Extractors = {}

-- Add an extractor for each output kind

Problem.Extractors[1] = {

    -- Export the mesh in a GMSH file

    kind = 'GMSH',
    writeAs = 'NodesElements',
    outputFile = 'pfem/fluid.msh',
    whatToWrite = {'p', 'velocity'},
    timeBetweenWriting = math.huge
}

Problem.Extractors[2] = {

    -- Export the total fluid mass

    kind = 'Global',
    whatToWrite = 'mass',
    outputFile = 'mass.txt',
    timeBetweenWriting = math.huge
}

Problem.Material = {}

-- First material is the fluid

Problem.Material[1] = {

    -- Parameters for the viscosity

    Stress = {
        type = 'NewtonianFluid',
        mu = 1e-3
    },
    
    -- Parameters for the fluid bulk

    StateEquation = {
        type = 'TaitMurnaghan',
        rho0 = 1000,
        K0 = 2.2e+7,
        K0p = 7.6,
        p0 = 0
    },

    -- Parameters for surface tension

    SurfaceStress = {
        type = 'SurfaceTension',
        gamma = 7e-2
    }
}

Problem.Solver = {

    -- Initial conditions and type

    IC = {},
    type = 'Explicit',
    timeIntegration = 'CDS',

    -- Factors of time step changes

    securityCoeff = 0.1,

    -- Enable or disable algorithms

    adaptDT = true,
    maxDT = math.huge,
    initialDT = math.huge
}

Problem.Solver.MomEq = {

    -- Enable the fluid-structure interface

    BC = {FSIVExt = true},

    -- Other simulation parameters

    pExt = 0,
    bodyForce = {0, -9.81}
}

Problem.Solver.ContEq = {

    -- Enable the boundary conditions

    BC = {},

    -- Define the solver algorithms

    version = 'DpDt',
    stabilization = 'CLS'
}

-- Initial Conditions

function Problem.Solver.IC.initStates(x, y, z)
	return {0, 0, 0, Problem.Material[1].StateEquation.rho0, 0, 0}
end

-- Momentum Equation BC

function Problem.Solver.MomEq.BC.ReservoirV(x, y, z, t)
	return 0, 0
end


