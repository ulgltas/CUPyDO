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
    
    mshFile = runPath+os.sep+'birdStrike_lsDyna_benchmark_Mtf_Pfem.msh'
    
    rho0 = 781.
    mu = 0.01
    U0 = 0.05
    V0 = -0.01
    N = 10

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
    print msh
    
    scheme = w.BackwardEuler(msh, pbl)

    w.Medium(msh, 15, mu, rho0, 3)
    w.Medium(msh, 17, mu, rho0, 1)
    
    # boundaries
    w.Boundary(msh, 14, 3, 0.0)
    w.Boundary(msh, 15, 1, 0.0)
    w.Boundary(msh, 15, 2, 0.0)
    
    #Initial velocity
    bird = w.Group(msh, 17)
    loadingset = w.LoadingSet(msh)
    loadingset.add(1,w.InitialVelocity(msh,bird,U0,V0,0.))
    
    scheme.savefreq=1
    scheme.nthreads=1
    scheme.gamma = 0.5
    scheme.omega = 0.5
    scheme.addRemoveNodesOption = True
    scheme.tollNLalgo = 1e-7
    
    #Results
    extManager = w.ExtractorsManager(msh)
    extManager.add(1,w.PositionExtractor(msh,4))
    extManager.add(2,w.PositionExtractor(msh,3))
    '''extManager.add(1,w.IntForceExtractor(msh,"Wall"))
    extManager.add(2,w.ExtForceExtractor(msh,"Wall"))
    extManager.add(3,w.IneForceExtractor(msh,"Wall"))
    extManager.add(4,w.PressureExtractor(msh,"Wall"))
    extManager.add(5,w.PositionExtractor(msh,"Wall"))
    extManager.add(6,wt.KineticEnergyExtractor(msh,pbl,"Body"))
    extManager.add(7,wt.ViscousEnergyExtractor(msh,pbl,scheme,"Body"))'''
    
    gui = v.MeshViewer(msh, scheme, True) 
    
    return Module(w, msh, pbl, bird, loadingset, scheme, extManager, gui)

def getRealTimeExtractorsList(pfem):
    
    extractorsList = []

    # --- Extractors list starts --- #
    # --- Extractors list ends --- #

    return extractorsList