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

#!/usr/bin/env python
# -*- coding: latin-1; -*-

# SU2Interface.py
# Python interface between the modal data base and CUPyDO.
# Authors D. THOMAS

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import numpy as np
from cupydo.genericSolvers import SolidSolver
from math import *

# ----------------------------------------------------------------------
#  Modal solver interface class
# ----------------------------------------------------------------------

class modalInterpreter(SolidSolver):
    """
    Description
    """

    def __init__(self, confFile, freq):
        """
        Des.
        """

        print("\n***************************** Initialize modal interpreter *****************************\n")

        self.modalDataFileName = confFile
        self.frequency = freq
        self.interfaceID = None
        self.nNodes = self.__linecount(self.modalDataFileName)-1
        self.nHaloNode = 0
        self.nPhysicalNodes = self.nNodes

        self.nodalMod_X = np.zeros(self.nPhysicalNodes, dtype=float)
        self.nodalMod_Y = np.zeros(self.nPhysicalNodes, dtype=float)
        self.nodalMod_Z = np.zeros(self.nPhysicalNodes, dtype=float)

        self.nodalDisp_X = np.zeros(self.nPhysicalNodes, dtype=float)
        self.nodalDisp_Y = np.zeros(self.nPhysicalNodes, dtype=float)
        self.nodalDisp_Z = np.zeros(self.nPhysicalNodes, dtype=float)

        self.nodalCoord_X = np.zeros(self.nPhysicalNodes, dtype=float)
        self.nodalCoord_Y = np.zeros(self.nPhysicalNodes, dtype=float)
        self.nodalCoord_Z = np.zeros(self.nPhysicalNodes, dtype=float)

        self.nodalGlobalIndex = np.zeros(self.nPhysicalNodes, dtype=int)

        self.haloNodeList = {}

        self.__readModalData()

    def __readModalData(self):
        """
        Des.
        """

        print('Counted modal data for {} nodes.'.format(self.nPhysicalNodes))
        iVertex = 0

        with open(self.modalDataFileName, 'r') as modalDataFile:
            print('Opened modal data file ' + self.modalDataFileName + '.')
            line = modalDataFile.readline()
            if line:
                pos = line.find('Global_Index')
                if pos == -1:
                    raise Exception("ERROR")
                while 1:
                    line = modalDataFile.readline()
                    if not line:
                        break
                    line = line.split(',')
                    for ii in range(len(line)):
                        line[ii] = line[ii].strip(' ')
                        line[ii] = line[ii].strip('\n\r')
                    self.nodalGlobalIndex[iVertex] = int(line[0])
                    self.nodalCoord_X[iVertex] = float(line[1])
                    self.nodalCoord_Y[iVertex] = float(line[2])
                    self.nodalCoord_Z[iVertex] = float(line[3])
                    self.nodalMod_X[iVertex] = float(line[4])
                    self.nodalMod_Y[iVertex] = float(line[5])
                    self.nodalMod_Z[iVertex] = float(line[6])
                    iVertex += 1

        if iVertex != self.nPhysicalNodes:
            raise("ERROR")
        else:
            print('Read modal data for {} nodes.'.format(iVertex))

    def __linecount(self, fileName):
        """
        Des
        """

        count = 0

        with open(fileName, 'r') as thefile:
            while 1:
                line = thefile.readline()
                if not line:
                    break
                count += 1

        return count

    def run(self, t1, t2):
        """
        Des.
        """

        mult = sin(2*3.141592*self.frequency*t2)

        self.nodalDisp_X = self.nodalMod_X*mult
        self.nodalDisp_Y = self.nodalMod_Y*mult
        self.nodalDisp_Z = self.nodalMod_Z*mult

    def getNodalInitialPositions(self):
        """
        Des.
        """

        return (self.nodalCoord_X, self.nodalCoord_Y, self.nodalCoord_Z)

    def getNodalIndex(self, iVertex):
        """
        Des.
        """

        if iVertex >= self.nPhysicalNodes or iVertex < 0:
            raise Exception('INDEX OUT OF BOUND.')
        else:
            return self.nodalGlobalIndex[iVertex]

    def exit(self):
        """
        Des.
        """

        print("***************************** Exit modal interpreter *****************************")


