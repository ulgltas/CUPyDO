#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# original name:

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


class Module(object):
    def __init__(self, w, msh, pbl, contactTag, solScheme,
                 nonLinAlgo, convCriterion, scheme,
                 extManager, bndno):
        self.w = w
        self.msh = msh
        self.pbl = pbl
        self.contactTag = contactTag
        self.solScheme = solScheme
        self.nonLinAlgo = nonLinAlgo
        self.convCriterion = convCriterion
        self.scheme = scheme
        self.extManager = extManager
        self.bndno = bndno


def getPfem():
    global w
    if w:
        return w
    w = pfem

    mshFile = runPath+os.sep+'matching.msh'

    rho0 = 1000.
    mu = 0.001

    pbl = w.Problem()
    pbl.rho0 = rho0
    pbl.mu = mu
    pbl.alpha = 1.2
    pbl.extP = 0.
    pbl.scalingU = 5.0
    pbl.bodyForceY = -9.81

    msh = w.MshData(pbl)
    msh.load(mshFile)
    print(msh)

    toll = 1e-6
    nItMax = 20

    solScheme = w.SchemeMonolithicPSPG(msh, pbl)
    convCriterion = w.ForceBalanceCriterion(msh, pbl, toll, 9.81*rho0)
    nonLinAlgo = w.PicardAlgorithm(solScheme, convCriterion, nItMax)

    scheme = w.TimeIntegration(msh, pbl, solScheme)

    bndno = "FSInterface"

    contactTag = w.Tag(100, "Contact", 2)
    msh.ptags[100] = contactTag
    msh.ntags["Contact"] = contactTag

    w.Medium(msh, "Contact", w.CONTACT, 0., 0.)
    w.Medium(msh, "FSInterface", w.SOLID, 0., 0.)
    w.Medium(msh, "Reservoir", w.MASTER_FLUID, mu, rho0)
    w.Medium(msh, "WaterColumn", w.MASTER_FLUID, mu, rho0)

    # boundaries
    w.Boundary(msh, "FreeSurface", 3, pbl.extP)
    w.Boundary(msh, "Reservoir", 1, 0.0)
    w.Boundary(msh, "Reservoir", 2, 0.0)
    w.Boundary(msh, "FSInterface", 1, 0.0)
    w.Boundary(msh, "FSInterface", 2, 0.0)

    scheme.savefreq = 100
    scheme.gamma = 0.6
    scheme.omega = 0.5
    scheme.addRemoveNodesOption = True

    # Results
    extManager = w.ExtractorsManager(msh)
    extManager.add(1, w.PositionExtractor(msh, 6))
    '''extManager.add(1,w.PositionExtractor(msh,1))
    extManager.add(2,w.IntForceExtractor(msh,1))
    extManager.add(3,w.ExtForceExtractor(msh,1))
    extManager.add(4,w.IneForceExtractor(msh,1))'''
    extManager.add(5, w.IntForceExtractor(msh, 6))
    '''extManager.add(6,w.ExtForceExtractor(msh,6))
    extManager.add(7,w.IneForceExtractor(msh,6))
    extManager.add(8,w.PressureExtractor(msh,1))
    extManager.add(9,w.VelocityExtractor(msh,"Water"))
    extManager.add(10,w.MassExtractor(msh,"Water"))
    extManager.add(11,wt.KineticEnergyExtractor(msh,pbl,"Water"))
    extManager.add(12,wt.ViscousEnergyExtractor(msh,pbl,scheme,"Water"))'''

    return Module(w, msh, pbl, contactTag, solScheme, nonLinAlgo, convCriterion, scheme, extManager, bndno)


def getRealTimeExtractorsList(pfem):
    return []
