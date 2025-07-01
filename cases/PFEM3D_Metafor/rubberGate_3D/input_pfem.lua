Problem = {

    mechanical = true,
    verboseOutput = true,
    autoRemeshing = false,
    simulationTime = math.huge
}

Problem.Mesh = {

    remeshAlgo = 'Tetgen_Edge',
    mshFile = 'geometryF.msh',
    deleteBoundElements = {'FSI'},
    localHcharGroups = {'Reservoir', 'FSI', 'FreeSurface'},
    exclusionZones = {{0.1, -0.05, 0.079, 0.2, 0.05, 0.14}},
    boundingBox = {0, -0.05, 0, 0.205, 0.05, 0.14},

    alpha = 1.7,
    omega = 0.4,
    gamma = 0.8,
    gammaFS = 0.8,
    hchar = 2.5e-3,
    omegaEdge = 1.6,
    epsilonLS = 0.15,
    minHeightFactor = 1e-2,

    addOnFS = true,
    useLevelSet = true,
    keepFluidElements = true,
    deleteFlyingNodes = false
}

Problem.Extractors = {}

Problem.Extractors[1] = {

    kind = 'GMSH',
    outputFile = 'pfem/output.msh',
    whatToWriteNode = {'p', 'velocity'},
    timeBetweenWriting = math.huge
}

Problem.Extractors[2] = {

    kind = 'Global',
    whatToWrite = 'mass',
    outputFile = 'mass.txt',
    timeBetweenWriting = math.huge
}

Problem.Material = {}

Problem.Material[1] = {

    Stress = {
        type = 'NewtonianFluid',
        mu = 1e-3
    },

    StateEquation = {
        type = 'Incompressible',
        rho = 1000
    },

    SurfaceStress = {
        type = 'SurfaceTension',
        gamma = 0
    }
}

Problem.Solver = {

    IC = {},
    type = 'Implicit',

    coeffDTDecrease = 2,
    coeffDTincrease = 1,

    adaptDT = true,
    maxDT = math.huge,
    initialDT = math.huge
}

Problem.Solver.MomContEq = {

    BC = {FSIVExt = true},

    nlAlgo = 'Picard',
    systemForm = 'FracStep',
    timeIntegration = 'BackwardEuler',
    residual = 'U_P',

    pExt = 0,
    maxIter = 25,
    gammaFS = 0.5,
    minRes = 1e-6,
    toleranceProj = 1e-16,
    toleranceVApp = 1e-16,
    toleranceVCorr = 1e-16,
    toleranceNeumann = 1e-16,
    bodyForce = {0, 0, -9.81}
}

function Problem.Solver.IC.initStates(x, y, z)
    return {0, 0, 0, 0}
end

function Problem.Solver.MomContEq.BC.ReservoirV(x, y, z, t)
    return 0, 0, 0
end

function Problem.Solver.MomContEq.BC.BottomV(x, y, z, t)
    return 0, 0, 0
end

function Problem.Mesh.computeHcharFromDistance(x, y, z, t, dist)

	local d = Problem.Mesh.hchar+math.pow(dist, 2)*4
    return math.min(d, 1e-2)
end
