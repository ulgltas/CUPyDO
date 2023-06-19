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
    def __init__(self, w, msh, pbl, bird, loadingset, scheme, extManager):
       self.w = w
       self.msh = msh
       self.pbl = pbl 
       self.bird = bird
       self.loadingset = loadingset      
       self.scheme = scheme
       self.extManager = extManager

def getPfem():
    global w
    if w: return w
    w = pfem
    
    mshFile = runPath+os.sep+'birdImpact_deformable_panel_Mtf_Pfem.msh'
    
    print('mshFile: ', mshFile)
    
    rho0 = 1000.
    mu = 0.001
    U0 = 100

    pbl = w.Problem()
    pbl.rho0 = rho0
    pbl.mu = mu
    pbl.nonLinAlgorithm = 1
    pbl.solScheme = 1
    pbl.alpha = 1.2
    pbl.beta = 0.
    pbl.gravity = 0.
    pbl.scalingU = U0

    msh = w.MshData(pbl)
    msh.load(mshFile)
    print(msh)
    
    scheme = w.BackwardEuler(msh, pbl)

    w.Medium(msh, 13, mu, rho0, 3)
    w.Medium(msh, 16, mu, rho0, 1)
    
    # boundaries
    w.Boundary(msh, 15, 3, 0.0)
    w.Boundary(msh, 13, 1, 0.0)
    w.Boundary(msh, 13, 2, 0.0)
    
    #Initial velocity
    bird = w.Group(msh, 16)
    loadingset = w.LoadingSet(msh)
    loadingset.add(1,w.InitialVelocity(msh,bird,0.,-U0,0.))

    scheme.gamma = 0.5
    scheme.omega = 0.5
    scheme.addRemoveNodesOption = True
    scheme.tollNLalgo = 1e-7
    
    #Results
    extManager = w.ExtractorsManager(msh)
    extManager.add(1,w.PositionExtractor(msh,310))
    extManager.add(6,wt.KineticEnergyExtractor(msh,pbl,16))
    extManager.add(7,wt.ViscousEnergyExtractor(msh,pbl,scheme,16))
    extManager.add(8,wt.PressureWorkExtractor(msh,pbl,scheme,16))
    extManager.add(9,w.MassExtractor(msh,pbl,16))

    return Module(w, msh, pbl, bird, loadingset, scheme, extManager)

def getRealTimeExtractorsList(pfem):
    
    extractorsList = []

    # --- Extractors list starts --- #
    # --- Extractors list ends --- #

    return extractorsList