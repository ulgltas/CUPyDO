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
    outputFile = 'pfem/fluid.msh',
    whatToWriteNode = {'p', 'velocity'},
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
    type = 'Implicit',

    -- Factors of time step changes

    coeffDTDecrease = 2,
    coeffDTincrease = 1,

    -- Enable or disable algorithms

    adaptDT = true,
    maxDT = math.huge,
    initialDT = math.huge
}

Problem.Solver.MomContEq = {

    -- Enable the fluid-structure interface

    BC = {FSIVExt = true},

    -- Define the solver algorithms

    nlAlgo = 'Picard',
    systemForm = 'Monolithic',
    timeIntegration = 'BackwardEuler',
    residual = 'Ax_f',

    -- Other simulation parameters

    pExt = 0,
    maxIter = 25,
    minRes = 1e-8,
    bodyForce = {0, -9.81}
}

-- Initial Conditions

function Problem.Solver.IC.initStates(x, y, z)
	return {0, 0, 0}
end
-- Momentum Continuity Equation BC

function Problem.Solver.MomContEq.BC.ReservoirV(x, y, z, t)
	return 0, 0
end


