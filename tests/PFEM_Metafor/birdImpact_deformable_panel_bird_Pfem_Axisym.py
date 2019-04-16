#! /usr/bin/env python
# -*- coding: latin-1; -*-

''' 

Copyright 2018 University of Li�ge

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 

'''

import sys, os, os.path

runPath = os.path.dirname(sys.modules[__name__].__file__)
filePath = os.path.abspath(os.path.dirname(__file__))
fileName = os.path.splitext(os.path.basename(__file__))[0]

import pfem

import pfemtools as wt
import viewer as v
    
w = None

class Module:
    def __init__(self, w, msh, pbl, solScheme, nonLinAlgo, convCriterion, bird, loadingset, scheme, extManager, gui):
       self.w = w
       self.msh = msh
       self.pbl = pbl       
       self.solScheme = solScheme
       self.nonLinAlgo = nonLinAlgo
       self.convCriterion = convCriterion
       self.bird = bird
       self.loadingset = loadingset
       self.scheme = scheme
       self.extManager = extManager
       self.gui = gui

def getPfem():
    global w
    if w: return w
    w = pfem
    
    mshFile = runPath+os.sep+'birdImpact_deformable_panel_Mtf_Pfem_Axisym.msh'
    
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
    pbl.isAxisymmetric = True
    
    msh = w.MshData(pbl)
    msh.load(mshFile)
    print msh
    
    toll = 1e-6
    nItMax = 10
    
    solScheme = w.SchemeMonolithicPSPG(msh, pbl)
    convCriterion = w.ForcesBalanceNormedBodyForceCriterion(msh, pbl, toll)
    nonLinAlgo = w.PicardAlgorithm(solScheme, convCriterion, nItMax)
    
    scheme = w.BackwardEuler(msh, pbl, nonLinAlgo)

    w.Medium(msh, 13, 0., 0., 3)
    w.Medium(msh, 17, mu, rho0, 1)
    
    # boundaries
    w.Boundary(msh, 15, 7, 0.0)
    w.Boundary(msh, 16, 3, 0.0)
    w.Boundary(msh, 13, 1, 0.0)
    w.Boundary(msh, 13, 2, 0.0)
    
    # --- Necessary for the correct elimination of nodes along the axis of symmetry ---
    n4 = msh.getNode(4)
    n4.isOnAxisOfSymmetry = True
    # ---
    
    #Initial velocity
    bird = w.Group(msh, 17)
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
    extManager.add(1,w.PositionExtractor(msh,4))
    extManager.add(6,wt.KineticEnergyExtractor(msh,pbl,17))
    extManager.add(7,wt.ViscousEnergyExtractor(msh,pbl,scheme,17))
    extManager.add(8,wt.PressureWorkExtractor(msh,pbl,scheme,17))
    extManager.add(9,w.MassExtractor(msh,pbl,17))
    
    gui = v.MeshViewer(msh, scheme, True) 
    
    return Module(w, msh, pbl, solScheme, nonLinAlgo, convCriterion, bird, loadingset, scheme, extManager, gui)

def getRealTimeExtractorsList(pfem):
    
    extractorsList = []

    # --- Extractors list starts --- #
    # --- Extractors list ends --- #

    return extractorsList