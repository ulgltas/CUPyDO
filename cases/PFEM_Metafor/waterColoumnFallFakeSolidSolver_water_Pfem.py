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
    def __init__(self, w, msh, pbl, scheme, extManager):
       self.w = w
       self.msh = msh
       self.pbl = pbl       
       self.scheme = scheme
       self.extManager = extManager

def getPfem():
    global w
    if w: return w
    w = pfem
    
    mshFile = runPath+os.sep+'waterColoumnFallFakeSolidSolver.msh'
    
    rho0 = 1000.
    mu = 0.001
    
    pbl = w.Problem()
    pbl.rho0 = rho0
    pbl.mu = mu
    pbl.nonLinAlgorithm = 1
    pbl.solScheme = 1
    pbl.alpha = 1.2
    pbl.extP = 100.
    pbl.scalingU = 1.0
    
    msh = w.MshData(pbl)
    msh.load(mshFile)
    print(msh)
    
    scheme = w.BackwardEuler(msh, pbl)
    
    w.Medium(msh, "Walls", mu, rho0, 3)
    w.Medium(msh, "Face1", mu, rho0, 1)
    
    # boundaries
    w.Boundary(msh, "Free", 3, pbl.extP)
    w.Boundary(msh, "Walls", 1, 0.0)
    w.Boundary(msh, "Walls", 2, 0.0)
    
    scheme.gamma = 0.5
    scheme.omega = 0.5
    scheme.addRemoveNodesOption = True
    scheme.tollNLalgo = 1e-12
    
    #Results
    extManager = w.ExtractorsManager(msh)
    extManager.add(1,w.PositionExtractor(msh,1))
    extManager.add(2,w.IntForceExtractor(msh,1))
    extManager.add(3,w.ExtForceExtractor(msh,1))
    extManager.add(4,w.IneForceExtractor(msh,1))
    extManager.add(5,w.IntForceExtractor(msh,6))
    extManager.add(6,w.ExtForceExtractor(msh,6))
    extManager.add(7,w.IneForceExtractor(msh,6))
    extManager.add(8,w.PressureExtractor(msh,1))
    extManager.add(9,w.VelocityExtractor(msh,"Face1"))
    extManager.add(10,w.MassExtractor(msh,pbl,"Face1"))
    extManager.add(11,wt.KineticEnergyExtractor(msh,pbl,"Face1"))
    extManager.add(12,wt.ViscousEnergyExtractor(msh,pbl,scheme,"Face1"))
    
    return Module(w, msh, pbl, scheme, extManager)