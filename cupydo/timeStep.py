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

Metafor.py
Python interface between the wrapper of Metafor and CUPyDO.
Authors R. BOMAN, M.L. CERQUAGLIA, D. THOMAS

'''

import math

# %% Coupling Time Step Manager

class TimeStep(object):
    def __init__(self,dt,dtSave):

        self.time = 0
        self.minDt = 1e-9
        self.division = int(2)
        self.maxDt = self.dt = dt
        self.next = self.dtSave = dtSave
        self.timeIter = 0

    def timeFrame(self):
        return self.time,self.time+self.dt

    # Update and return if save is required

    def mustSave(self):

        if self.dtSave <= 0: return True

        output = (self.time >= self.next)
        next = math.floor(self.time/self.dtSave)
        self.next = (next+1)*self.dtSave
        return output
        
    # Update the coupling time step

    def update(self,verified):

        if not verified:
            
            self.dt /= self.division
            if self.dt < self.minDt:
                raise Exception('Reached minimal time step')

        else:

            self.timeIter += 1
            self.time += self.dt
            self.dt = math.pow(self.division,1/7)*self.dt
            self.dt = min(self.dt,self.maxDt)
