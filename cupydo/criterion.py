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
from .interfaceData import FlexInterfaceData

# ----------------------------------------------------------------------
#    Criterion class
# ----------------------------------------------------------------------

class Criterion(object):

    def __init__(self, p):

        self.epsilon = 0
        self.epsilonCHT = 0
        self.tol = p['tol']

    # Reset all the convergence indicators

    def reset(self):

        self.epsilon = 0
        self.epsilonCHT = 0

    def isVerified(self):
        return (self.epsilon < self.tol) and (self.epsilonCHT < self.tol)

class NormCriterion(Criterion):
    
    def __init__(self, p):
        Criterion.__init__(self, p)

    # Update the mechanical residual

    def update(self, prediction, residual):

        res = np.array(residual.norm())
        res /= (np.array(prediction.norm())+self.tol)
        self.epsilon = np.linalg.norm(res)

    # Update ththermal residual

    def update_CHT(self, prediction, residual):

        res = np.array(residual.norm())
        res /= (np.array(prediction.norm())+self.tol)
        self.epsilonCHT = np.linalg.norm(res)
