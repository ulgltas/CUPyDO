Problem = {

    -- Main simulation parameters

    thermal = true,
    mechanical = true,
    boussinesq = true,
    verboseOutput = true,
    autoRemeshing = false,
    simulationTime = math.huge
}

Problem.Mesh = {

    -- Input mesh and bounding box

    remeshAlgo = 'CGAL',
    mshFile = 'geometryF.msh',
    deleteBoundElements = {'FSI'},
    boundingBox = {0, 0, 0.9, 1},
    exclusionZones = {},

    -- Remeshing internal parameters

    alpha = 1.2,
    omega = 0.5,
    gamma = 0.6,
    hchar = 0.25/33,
    gammaFS = 0.2,
    epsilonLS = 0.1,
    gammaEdge = 0.2,
    minHeightFactor = 1e-3,

    -- Enable or disable algorithms

    addOnFS = false,
    useLevelSet = true,
    keepFluidElements = true,
    deleteFlyingNodes = true
}

Problem.Extractors = {}

-- Add an extractor for each output kind

Problem.Extractors[1] = {

    -- Export the mesh in a GMSH file

    kind = 'GMSH',
    outputFile = 'pfem/fluid.msh',
    whatToWriteNode = {'T', 'velocity'},
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
        mu = 5e-3
    },
    
    -- Parameters for the fluid bulk

    StateEquation = {
        type = 'Incompressible',
        rho = 1000
    },

    -- Parameters for surface tension

    SurfaceStress = {
        type = 'SurfaceTension',
        gamma = 0
    },

    -- Parameters for heat capacity

    CaloricStateEq = {
        type = 'LinearHeatCapacity',
        DcpDT = 0,
        Tr = 340,
        cp = 1e3
    },

    -- Parameters for Fourier flux

    HeatFlux = {
        type = 'LinearFourierFlux',
        DkDT = 0,
        Tr = 340,
        k = 20
    },

    -- Parameters for cooling law
    
    CoolingLaw = {
        type = 'LinearCoolingLaw',
        Tinf = 340,
        h = 0
    },

    -- Parameters for optical material
    
    OpticalProperties = {
        type = 'ConstantOpticalProperties',
        absorptivity = 0,
        emissivity = 0,
        sigmaSB = 0,
        Tinf = 340
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
    initialDT = math.huge,
    solveHeatFirst = false
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

Problem.Solver.HeatEq = {

    -- Enable the fluid-structure interface

    BC = {FSITExt = true},

    -- Define the solver algorithms

    nlAlgo = 'Picard',
    timeIntegration = 'BackwardEuler',
    residual = 'Ax_f',

    -- Other simulation parameters

    maxIter = 25,
    minRes = 1e-6,
    tolerance = 1e-16,
}

-- Initial Conditions

function Problem.Solver.IC.initStates(x, y, z)
	return {0, 0, 0, 340}
end

-- Momentum Continuity Equation BC

function Problem.Solver.MomContEq.BC.WallV(x, y, z, t)
	return 0, 0
end

-- Heat Equation BC

function Problem.Solver.HeatEq.BC.WallQ(x, y, z, t)
    return 0, 0
end

function Problem.Solver.HeatEq.BC.FreeSurfaceQ(x, y, z, t)
    return 0, 0
end
