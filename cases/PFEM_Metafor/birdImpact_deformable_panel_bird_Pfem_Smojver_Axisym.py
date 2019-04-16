#! /usr/bin/env python
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
    def __init__(self, w, msh, pbl, bird, loadingset, scheme, extManager, gui):
       self.w = w
       self.msh = msh
       self.pbl = pbl
       self.bird = bird
       self.loadingset = loadingset       
       self.scheme = scheme
       self.extManager = extManager
       self.gui = gui

def getPfem():
    global w
    if w: return w
    w = pfem
    
    mshFile = runPath+os.sep+'birdImpact_deformable_panel_Mtf_Pfem_Smojver_Axisym.msh'
    
    rho0 = 938.
    mu = 0.001
    U0 = 145.7

    pbl = w.Problem()
    pbl.rho0 = rho0
    pbl.mu = mu
    pbl.nonLinAlgorithm = 1
    pbl.solScheme = 1
    pbl.alpha = 1.2
    pbl.beta = 0.
    pbl.gravity = 0.
    pbl.scalingU = U0
    pbl.isAxisymmetric = True
    
    msh = w.MshData(pbl)
    msh.load(mshFile)
    print msh
    
    scheme = w.BackwardEuler(msh, pbl)

    w.Medium(msh, 17, mu, rho0, 3)
    w.Medium(msh, 21, mu, rho0, 1)
    
    # boundaries
    w.Boundary(msh, 20, 7, 0.0)
    w.Boundary(msh, 19, 3, 0.0)
    w.Boundary(msh, 17, 1, 0.0)
    w.Boundary(msh, 17, 2, 0.0)
    
    # --- Necessary for the correct elimination of nodes along the axis of symmetry ---
    n4 = msh.getNode(6)
    n4.isOnAxisOfSymmetry = True
    # ---
    
    #Initial velocity
    bird = w.Group(msh, 21)
    loadingset = w.LoadingSet(msh)
    loadingset.add(1,w.InitialVelocity(msh,bird,0.,-U0,0.))
    
    scheme.savefreq=1
    scheme.nthreads=3
    scheme.gamma = 0.5
    scheme.omega = 0.5
    scheme.addRemoveNodesOption = True
    scheme.tollNLalgo = 1e-7
    
    #Results
    extManager = w.ExtractorsManager(msh)
    extManager.add(1,w.PositionExtractor(msh,6))
    extManager.add(6,wt.KineticEnergyExtractor(msh,pbl,21))
    extManager.add(7,wt.ViscousEnergyExtractor(msh,pbl,scheme,21))
    extManager.add(8,wt.PressureWorkExtractor(msh,pbl,scheme,21))
    extManager.add(9,w.MassExtractor(msh,pbl,21))
    
    gui = v.MeshViewer(msh, scheme, True) 
    
    return Module(w, msh, pbl, bird, loadingset, scheme, extManager, gui)

def getRealTimeExtractorsList(pfem):
    
    extractorsList = []

    # --- Extractors list starts --- #
    # --- Extractors list ends --- #

    return extractorsList