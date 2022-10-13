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
    def __init__(self, w, msh, pbl, solScheme, nonLinAlgo, convCriterion, scheme, extManager, gui):
       self.w = w
       self.msh = msh
       self.pbl = pbl       
       self.solScheme = solScheme
       self.nonLinAlgo = nonLinAlgo
       self.convCriterion = convCriterion
       self.scheme = scheme
       self.extManager = extManager
       self.gui = gui

def getPfem():
    global w
    if w: return w
    w = pfem
    
    mshFile = runPath+os.sep+'waterColumnOnElasticColumn_Mtf_Pfem.msh'
    
    rho0 = 1000.
    mu = 0.
    
    pbl = w.Problem()
    pbl.rho0 = rho0
    pbl.mu = mu
    pbl.alpha = 2.0
    pbl.extP = 0.
    pbl.scalingU = 15.0
    pbl.bodyForceY = -10.
    pbl.remeshing = False
    
    msh = w.MshData(pbl)
    msh.load(mshFile)
    print msh
    
    toll = 1e-6
    nItMax = 20
    
    solScheme = w.SchemeMonolithicPSPG(msh, pbl)
    solScheme.tauPSPGcomputation = 3
    convCriterion = w.ForcesBalanceNormedBodyForceCriterion(msh, pbl, toll)
    # convCriterion = w.PositionIncrementCriterion(msh, pbl, toll)
    nonLinAlgo = w.PicardAlgorithm(solScheme, convCriterion, nItMax)
    
    scheme = w.BackwardEuler(msh, pbl, nonLinAlgo)
    
    w.Medium(msh, 14, mu, rho0, 3)
    w.Medium(msh, 18, mu, rho0, 1)
    
    # boundaries
    w.Boundary(msh, 16, 3, pbl.extP)
    w.Boundary(msh, 15, 1, 0.0)
    w.Boundary(msh, 14, 1, 0.0)
    w.Boundary(msh, 14, 2, 0.0)
    
    scheme.nthreads = 1
    
    #Results
    extManager = w.ExtractorsManager(msh)
    extManager.add(1,w.PositionExtractor(msh,3))
    extManager.add(2,w.PositionExtractor(msh,6))
    extManager.add(3,w.PositionExtractor(msh,76))
    
    gui = v.MeshViewer(msh, scheme, True) 
    
    return Module(w, msh, pbl, solScheme, nonLinAlgo, convCriterion, scheme, extManager, gui)

def getRealTimeExtractorsList(pfem):
    
    extractorsList = []

    # --- Extractors list starts --- #
    # --- Extractors list ends --- #

    return extractorsList