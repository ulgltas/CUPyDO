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
    def __init__(self, w, msh, pbl, contactTag, solScheme, nonLinAlgo, convCriterion, scheme, extManager, gui):
       self.w = w
       self.msh = msh
       self.pbl = pbl       
       self.contactTag = contactTag
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
    
    mshFile = runPath+os.sep+'elasticContainerFilling.msh'
    
    rho0 = 1000.
    mu = 100.0
    
    pbl = w.Problem()
    pbl.rho0 = rho0
    pbl.mu = mu
    pbl.alpha = 1.2
    pbl.extP = 0.
    pbl.scalingU = 10.0
    pbl.bodyForceY = -9.81
    
    msh = w.MshData(pbl)
    msh.load(mshFile)
    print msh
    
    toll = 1e-6
    nItMax = 20
    
    solScheme = w.SchemeMonolithicPSPG(msh, pbl)
    convCriterion = w.ForcesBalanceNormedBodyForceCriterion(msh, pbl, toll)
    nonLinAlgo = w.PicardAlgorithm(solScheme, convCriterion, nItMax)
    
    scheme = w.BackwardEuler(msh, pbl, nonLinAlgo)
    
    contactTag = w.Tag(100, "Contact" , 2)
    msh.ptags[100] = contactTag
    msh.ntags["Contact"] = contactTag
    
    w.Medium(msh, 100, 0., 0., 0., 4)
    w.Medium(msh, 2, 0., 0., 3)
    w.Medium(msh, 4, 0., 0., 3)
    w.Medium(msh, 5, 0., 0., 3)
    w.Medium(msh, 6, mu, rho0, 1)
    
    # boundaries
    w.Boundary(msh, 1, 3, pbl.extP)
    w.Boundary(msh, 2, 1, 0.0)
    w.Boundary(msh, 2, 2, 0.0)
    w.Boundary(msh, 4, 1, 0.0)
    w.Boundary(msh, 4, 2, 0.0)
    w.Boundary(msh, 5, 1, 0.0)
    w.Boundary(msh, 5, 2, 0.0)
    
    scheme.savefreq=1
    scheme.nthreads=4
    scheme.gamma = 0.6
    scheme.omega = 0.5
    scheme.addRemoveNodesOption = True
    
    #Results
    extManager = w.ExtractorsManager(msh)
    extManager.add(1,w.PositionExtractor(msh,462))
    extManager.add(2,w.PositionExtractor(msh,463))
    
    gui = v.MeshViewer(msh, scheme, True) 
    
    return Module(w, msh, pbl, contactTag, solScheme, nonLinAlgo, convCriterion, scheme, extManager, gui)
    
def getRealTimeExtractorsList(pfem):
    
    extractorsList = []

    # --- Extractors list starts --- #
    # --- Extractors list ends --- #

    return extractorsList