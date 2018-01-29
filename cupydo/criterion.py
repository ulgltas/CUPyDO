#! /usr/bin/env python
# -*- coding: latin-1; -*-

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

criterion.py
Defines coupling criteria for CUPyDO.
Authors : David THOMAS, Marco Lucio CERQUAGLIA, Romain BOMAN

'''

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

from math import *

# ----------------------------------------------------------------------
#    Criterion class
# ----------------------------------------------------------------------

class Criterion:
    """
    Description
    """

    def __init__(self, tolerance, heatFluxTolerance = 1e12):
        """
        Description.
        """

        self.tol = tolerance
        self.epsilon = 0

        self.tolHeatFlux = heatFluxTolerance
        self.epsilonHeatFlux = 0

    def isVerified(self, epsilon, epsilonHeatFlux=0):
        """
        Description.
        """

        verifList = [False, False]

        if epsilon < self.tol:
            verifList[0] = True

        if epsilonHeatFlux < self.tolHeatFlux:
            verifList[1] = True

        if False in verifList:
            return False
        else:
            return True

class DispNormCriterion(Criterion):
    """
    Description.
    """

    def __init__(self, tolerance, heatFluxTolerance = 1e12):
        """
        Description.
        """

        Criterion.__init__(self, tolerance, heatFluxTolerance)

    def update(self, res):
        """
        Des.
        """

        normX, normY, normZ = res.norm()

        norm = sqrt(normX**2 + normY**2 + normZ**2)

        self.epsilon = norm

        return self.epsilon

    def updateHeatFlux(self, resHeatFlux):
        """
        Des.
        """

        if resHeatFlux != None:
            #normX, normY, normZ = resHeatFlux.norm()
            normList = resHeatFlux.norm()
            normSquare = 0.0
            for ii in range(len(normList)):
                normSquare += normList[ii]**2

            #norm = sqrt(normX**2 + normY**2 + normZ**2)
            norm = sqrt(normSquare)
        else:
            norm = 1.0

        self.epsilonHeatFlux = norm

        return self.epsilonHeatFlux
