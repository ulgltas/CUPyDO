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

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import math

# ----------------------------------------------------------------------
#    Time Step Manager Class
# ----------------------------------------------------------------------

class TimeStep(object):
    def __init__(self,Manager, FluidSolver, SolidSolver, p, mpiComm):

        self.mpiComm = mpiComm

        if self.mpiComm != None:
            self.myid = self.mpiComm.Get_rank()
            self.mpiSize = self.mpiComm.Get_size()
        else:
            self.myid = 0
            self.mpiSize = 1

        self.manager = Manager
        self.FluidSolver = FluidSolver
        self.SolidSolver = SolidSolver

        # Time step internal variables

        self.time = 0
        self.timeIter = 0
        self.minDt = 1e-9
        self.division = int(2)
        self.maxDt = self.dt = p['dt']
        self.next = self.dtSave = p['dtSave']

    def timeFrame(self):
        """
        Get [t1, t2] for fluid and solid run
        """
        return self.time,self.time+self.dt

    def updateSave(self):
        """
        Update next save time and export results if needed
        """

        if self.time >= self.next:
            self.FluidSolver.save(self.timeIter)
            if self.myid in self.manager.getSolidSolverProcessors():
                self.SolidSolver.save()

        if self.dtSave > 0:
            next = math.floor(self.time/self.dtSave)
            self.next = (next+1)*self.dtSave

    def updateTime(self,verified):
        """
        Update the current coupling time step
        """

        if not verified:
            
            self.dt /= self.division
            if self.dt < self.minDt:
                raise Exception('Reached minimal time step')

        else:
            
            self.timeIter += 1
            self.time += self.dt
            self.dt = math.pow(self.division,1/7)*self.dt
            self.dt = min(self.dt,self.maxDt)
