#! /usr/bin/env python3
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
    def __init__(self, w, msh, pbl, solScheme, nonLinAlgo, convCriterion, scheme, extManager):
       self.w = w
       self.msh = msh
       self.pbl = pbl
       self.solScheme = solScheme
       self.nonLinAlgo = nonLinAlgo
       self.convCriterion = convCriterion
       self.scheme = scheme
       self.extManager = extManager

def getPfem():
    global w
    if w: return w
    w = pfem
    
    mshFile = runPath+os.sep+'StaticCylinder.msh'
    
    rho0 = 1.0
    mu = 0.01
    
    pbl = w.Problem()
    pbl.rho0 = rho0
    pbl.mu = mu
    pbl.alpha = 1.4
    pbl.extP = 1.0
    pbl.scalingU = 1.0
    pbl.bodyForceY = -9.81
    
    msh = w.MshData(pbl)
    msh.load(mshFile)
    print(msh)
    
    toll = 1e-6
    nItMax = 20
    
    solScheme = w.SchemeMonolithicPSPG(msh, pbl)
    convCriterion = w.ForcesBalanceNormedBodyForceCriterion(msh, pbl, toll)
    nonLinAlgo = w.PicardAlgorithm(solScheme, convCriterion, nItMax)
    
    scheme = w.BackwardEuler(msh, pbl, nonLinAlgo)
    
    msh.ptags[9].name = "Cylinder"
    msh.ntags["Cylinder"] = msh.ptags[9]
    
    w.Medium(msh, 9, 0., 0., 3)
    w.Medium(msh, 14, mu, rho0, 1)
    
    # boundaries
    w.Boundary(msh, 8, 3, pbl.extP)
    w.Boundary(msh, 7, 1, 0.0)
    w.Boundary(msh, 7, 2, 0.0)
    w.Boundary(msh, 9, 1, 0.0)
    w.Boundary(msh, 9, 2, 0.0)

    scheme.gamma = 0.6
    scheme.omega = 0.5
    scheme.addRemoveNodesOption = True
    
    #Results
    extManager = w.ExtractorsManager(msh)
    extManager.add(1,w.PositionExtractor(msh,6))
    extManager.add(20,w.IntForceExtractor(msh,"Cylinder"))
    extManager.add(21,wt.TotalIntForceExtractor(msh,"Cylinder"))

    
    return Module(w, msh, pbl, solScheme, nonLinAlgo, convCriterion, scheme, extManager)
    
def getRealTimeExtractorsList(pfem):
    
    extractorsList = []

    # --- Extractors list starts --- #
    # --- Extractors list ends --- #

    return extractorsList