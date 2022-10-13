#! /usr/bin/env python3
# -*- coding: utf8 -*-

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

class Criterion(object):
    """
    Description
    """

    def __init__(self, tolerance, thermalTolerance = 1e12):
        """
        Description.
        """

        self.tol = tolerance
        self.epsilon = 0

        self.tolthermal = thermalTolerance
        self.epsilonThermal = 0

    def isVerified(self, epsilon, epsilonThermal=0):
        """
        Description.
        """

        verifList = [False, False]

        if epsilon < self.tol:
            verifList[0] = True

        if epsilonThermal < self.tolthermal:
            verifList[1] = True

        if False in verifList:
            return False
        else:
            return True

class DispNormCriterion(Criterion):
    """
    Description.
    """

    def __init__(self, tolerance, thermalTolerance = 1e12):
        """
        Description.
        """

        Criterion.__init__(self, tolerance, thermalTolerance)

    def update(self, res):
        """
        Des.
        """

        normX, normY, normZ = res.norm()

        norm = sqrt(normX**2 + normY**2 + normZ**2)

        self.epsilon = norm

        return self.epsilon

    def updateThermal(self, resThermal):
        """
        Des.
        """

        if resThermal != None:
            #normX, normY, normZ = resHeatFlux.norm()
            normList = resThermal.norm()
            normSquare = 0.0
            for ii in range(len(normList)):
                normSquare += normList[ii]**2

            #norm = sqrt(normX**2 + normY**2 + normZ**2)
            norm = sqrt(normSquare)
        else:
            norm = 1.0

        self.epsilonThermal = norm

        return self.epsilonThermal
