#! /usr/bin/env python
# -*- coding: utf-8 -*-
# original name: birdStrike_lsDyna_benchmark_bird_Pfem

''' 

Copyright 2018 University of Liège

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

import sys
import os

runPath = os.path.dirname(sys.modules[__name__].__file__)
filePath = os.path.abspath(os.path.dirname(__file__))
fileName = os.path.splitext(os.path.basename(__file__))[0]

import pfem
import pfem.tools.pfemtools as wt

w = None


class Module:
    def __init__(self, w, msh, pbl, solScheme, nonLinAlgo,
                 convCriterion, bird, loadingset, scheme, extManager, gui, bndno):
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
        self.bndno = bndno


def getPfem():
    global w
    if w:
        return w
    w = pfem

    mshFile = runPath + os.sep + 'lsdyna.msh'

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

    toll = 1e-6
    nItMax = 10

    solScheme = w.SchemeMonolithicPSPG(msh, pbl)
    convCriterion = w.ForcesBalanceNormedBodyForceCriterion(msh, pbl, toll)
    nonLinAlgo = w.PicardAlgorithm(solScheme, convCriterion, nItMax)

    scheme = w.BackwardEuler(msh, pbl, nonLinAlgo)

    bndno = 15 # fsi boundary

    w.Medium(msh, 15, 0., 0., 3)
    w.Medium(msh, 17, mu, rho0, 1)

    # boundaries
    w.Boundary(msh, 14, 3, 0.0)
    w.Boundary(msh, 15, 1, 0.0)
    w.Boundary(msh, 15, 2, 0.0)

    # Initial velocity
    bird = w.Group(msh, 17)
    loadingset = w.LoadingSet(msh)
    loadingset.add(1, w.InitialVelocity(msh, bird, U0, V0, 0.))

    scheme.savefreq = 1
    scheme.gamma = 0.5
    scheme.omega = 0.5
    scheme.addRemoveNodesOption = True
    scheme.tollNLalgo = 1e-7

    # Results
    extManager = w.ExtractorsManager(msh)
    extManager.add(1, w.PositionExtractor(msh, 4))
    extManager.add(2, w.PositionExtractor(msh, 3))
    '''extManager.add(1,w.IntForceExtractor(msh,"Wall"))
    extManager.add(2,w.ExtForceExtractor(msh,"Wall"))
    extManager.add(3,w.IneForceExtractor(msh,"Wall"))
    extManager.add(4,w.PressureExtractor(msh,"Wall"))
    extManager.add(5,w.PositionExtractor(msh,"Wall"))
    extManager.add(6,wt.KineticEnergyExtractor(msh,pbl,"Body"))
    extManager.add(7,wt.ViscousEnergyExtractor(msh,pbl,scheme,"Body"))'''

    import pfem.tools.link2vtk as v
    gui = v.Link2VTK(msh, scheme, True)

    return Module(w, msh, pbl, solScheme, nonLinAlgo,
                  convCriterion, bird, loadingset,
                  scheme, extManager, gui, bndno)


def getRealTimeExtractorsList(pfem):
    return []
