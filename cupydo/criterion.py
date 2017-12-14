#!/usr/bin/env python
# -*- coding: latin-1; -*-
#
# criterion.py
# Defines coupling criteria for CUPyDO.
# Authors : David THOMAS, Marco Lucio CERQUAGLIA, Romain BOMAN
#
# COPYRIGHT (C) University of Liège, 2017.

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
