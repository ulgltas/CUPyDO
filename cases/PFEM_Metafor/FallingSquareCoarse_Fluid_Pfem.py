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
    
    mshFile = runPath+os.sep+'MovingSquareCoarse.msh'
    
    # Physical parameters for the fluid
    rho0 = 1. #kg/m3
    mu = 0.
    
    # Problem definition
    pbl = w.Problem()
    pbl.rho0 = rho0
    pbl.mu = mu
    pbl.alpha = 1.4                         # Alpha-shapes parameter
    pbl.extP = 1.                           # External pressure (arbitrary in incompressible formulations)
    pbl.scalingU = 0.03                        # Global scaling velocity for PSPG stabilization
    pbl.bodyForceY = -10.                  # Gravity
    
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
    nonLinAlgo = w.PicardAlgorithm(solScheme, convCriterion, nItMax)
    
    # Time integration scheme
    scheme = w.BackwardEuler(msh, pbl, nonLinAlgo)
    
    msh.ptags[3].name = "Square"
    msh.ntags["Square"] = msh.ptags[3]
    
    w.Medium(msh, 3, 0., 0., 3)
    w.Medium(msh, 4, mu, rho0, 1)
    
    # boundaries
    w.Boundary(msh, 1, 1, 0.0)
    w.Boundary(msh, 1, 2, 0.0)
    w.Boundary(msh, 2, 3, pbl.extP)
    w.Boundary(msh, 3, 1, 0.0)
    w.Boundary(msh, 3, 2, 0.0)
    
    # Time integration
    scheme.nthreads=4
    scheme.gamma = 0.5
    scheme.omega = 0.5
    scheme.addRemoveNodesOption = True
    
    extManager = w.ExtractorsManager(msh)
    extManager.add(1,w.PositionExtractor(msh,"Square"))
    extManager.add(2,w.VelocityExtractor(msh,"Square"))
    extManager.add(3,w.IntForceExtractor(msh,"Square"))
    extManager.add(4,wt.TotalIntForceExtractor(msh,"Square"))
    extManager.add(5,wt.TotalExtForceExtractor(msh,"Square"))
    extManager.add(6,wt.TotalIneForceExtractor(msh,"Square"))
    extManager.add(7,w.PressureExtractor(msh,"Square"))
    
    return Module(w, msh, pbl, scheme, solScheme, nonLinAlgo, convCriterion, extManager)
    
def getRealTimeExtractorsList(pfem):
    
    extractorsList = []

    # --- Extractors list starts --- #
    # --- Extractors list ends --- #

    return extractorsList
