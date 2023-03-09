#! /usr/bin/env python3
# -*- coding: latin-1; -*-
# $Id: $

import sys, os, os.path

runPath = os.path.dirname(sys.modules[__name__].__file__)
filePath = os.path.abspath(os.path.dirname(__file__))
fileName = os.path.splitext(os.path.basename(__file__))[0]

import pfem

import pfemtools as wt
import viewer as v
    
w = None

class Module:
    def __init__(self, w, msh, pbl, scheme, solScheme, nonLinAlgo, convCriterion, extManager):
        self.w = w
        self.msh = msh
        self.pbl = pbl       
        self.scheme = scheme
        self.solScheme = solScheme
        self.nonLinAlgo = nonLinAlgo
        self.convCriterion = convCriterion
        self.extManager = extManager

def getPfem():
    global w
    if w: return w
    w = pfem
    
    mshFile = runPath+os.sep+'FallingCylinder_Franci.msh'
    
    # Physical parameters for the fluid
    rho0 = 1000. #kg/m3
    mu = 0.1
    
    # Problem definition
    pbl = w.Problem()
    pbl.rho0 = rho0
    pbl.mu = mu
    pbl.alpha = 2.0                         # Alpha-shapes parameter
    pbl.extP = 735.                           # External pressure (arbitrary in incompressible formulations)
    pbl.scalingU = 0.03                        # Global scaling velocity for PSPG stabilization
    pbl.bodyForceY = -9.81                  # Gravity
    
    # Initial mesh
    msh = w.MshData(pbl)
    msh.load(mshFile)
    print(msh)
    
    # Formulation (PSPG stabilization/Fractional step/Algebraic splitting)
    solScheme = w.SchemeMonolithicPSPG(msh, pbl)
    
    # Non-linear iteration algorithm
    toll = 1e-6
    nItMax = 20
    convCriterion = w.ForcesBalanceNormedBodyForceCriterion(msh, pbl, toll)
    # convCriterion = w.PositionIncrementCriterion(msh, pbl, toll)
    nonLinAlgo = w.PicardAlgorithm(solScheme, convCriterion, nItMax)
    
    # Time integration scheme
    scheme = w.BackwardEuler(msh, pbl, nonLinAlgo)
    
    w.Medium(msh, 15, 0., 0., 3)
    w.Medium(msh, 16, mu, rho0, 1)
    
    # boundaries
    w.Boundary(msh, 11, 2, 0.0)
    w.Boundary(msh, 12, 1, 0.0)
    w.Boundary(msh, 13, 2, 0.0)
    w.Boundary(msh, 14, 1, 0.0)
    w.Boundary(msh, 15, 1, 0.0)
    w.Boundary(msh, 15, 2, 0.0)
    
    # Time integration
    scheme.nthreads=4
    scheme.gamma = 0.5
    scheme.omega = 0.5
    scheme.addRemoveNodesOption = True
    
    extManager = w.ExtractorsManager(msh)
    extManager.add(1,w.PositionExtractor(msh,15))
    extManager.add(2,w.VelocityExtractor(msh,15))
    extManager.add(3,w.PressureExtractor(msh,15))
    extManager.add(4,w.IntForceExtractor(msh,15))
    
    return Module(w, msh, pbl, scheme, solScheme, nonLinAlgo, convCriterion, extManager)
    
def getRealTimeExtractorsList(pfem):
    
    extractorsList = []

    # --- Extractors list starts --- #
    # --- Extractors list ends --- #

    return extractorsList
