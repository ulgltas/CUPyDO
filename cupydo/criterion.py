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
import numpy as np

# ----------------------------------------------------------------------
#    Criterion class
# ----------------------------------------------------------------------

class Criterion(object):

    def __init__(self, p):

        self.mechanical = p['mechanical']
        self.thermal = p['thermal']
        self.reset()

        if self.mechanical:
            self.tol = p['mechanicalTol']

        if self.thermal:
            self.tolCHT = p['thermalTol']

    # Reset all the convergence indicators

    def reset(self):

        self.epsilon = np.inf
        self.epsilonCHT = np.inf

    def isVerified(self):

        verified = list()
        if self.mechanical:
            verified.append(self.epsilon < self.tol)
        if self.thermal:
            verified.append(self.epsilonCHT < self.tolCHT)

        return np.all(verified)

class RelativeCriterion(Criterion):
    
    def __init__(self, p):
        Criterion.__init__(self, p)

    # Update the mechanical residual

    def update(self, residual, prediction):

        res = np.array(residual.norm())
        res /= (np.array(prediction.norm())+self.tol)
        self.epsilon = np.linalg.norm(res)

    # Update the thermal residual

    def update_CHT(self, residual, prediction):

        res = np.array(residual.norm())
        res /= (np.array(prediction.norm())+self.tolCHT)
        self.epsilonCHT = np.linalg.norm(res)

class NormCriterion(Criterion):
    
    def __init__(self, p):
        Criterion.__init__(self, p)

    # Update the mechanical residual

    def update(self, residual, prediction):

        res = np.array(residual.norm())
        self.epsilon = np.linalg.norm(res)

    # Update the thermal residual

    def update_CHT(self, residual, prediction):

        res = np.array(residual.norm())
        self.epsilonCHT = np.linalg.norm(res)