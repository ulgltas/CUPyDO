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

GetDP.py
Python interface between the wrapper of GetDP and CUPyDO.
Authors C. GEUZAINE, D. THOMAS

'''

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import math
import os
import numpy as np
from ..genericSolvers import SolidSolver
from .getdp import *

# ----------------------------------------------------------------------
#  GetDP solver interface class
# ----------------------------------------------------------------------

class GetDP(SolidSolver):

    def __init__(self, p):


        self.testname = p['csdFile']
        self.pythonFlag = p['pythonFlag']
        self.regime = p['regime']
        self.resolution = p['resolution']
        self.currentDT = 1.0
        self.pathToGetDP = "/home/dthomas/InstalledSoftware/GetDP/bin/getdp"

        if self.pythonFlag:
            #from getdp import *
            GetDPSetNumber("Initialize", 1)
            GetDPSetNumber("OutputFiles", 1)
            GetDP(["getdp", self.testname, "-solve", self.resolution])
            GetDPPos = GetDPGetNumber("nodalPosition") #returns a std vector
            self.nNodes = int(GetDPPos[0])
            self.nHaloNode = 0
            self.nPhysicalNodes = int(GetDPPos[0])
            GetDPDomainDisp = GetDPGetNumber("nodalDisplacement")
            self.__nDomainNodes = int(GetDPDomainDisp[0])
            SolidSolver.__init__(self)
            self.nodalInterfIndex = self.__extractIndex(GetDPPos, 3)
            self.nodalDomainIndex = self.__extractIndex(GetDPDomainDisp, 3)
            self.nodalInitialPos_X = np.zeros(self.nPhysicalNodes)
            self.nodalInitialPos_Y = np.zeros(self.nPhysicalNodes)
            self.nodalInitialPos_Z = np.zeros(self.nPhysicalNodes)
            self.nodalInitialPos_X, self.nodalInitialPos_Y, self.nodalInitialPos_Z = self.__vecToVecArray(GetDPPos)
            self.__nodalDomainDisp_X = np.zeros(self.__nDomainNodes)
            self.__nodalDomainDisp_Y = np.zeros(self.__nDomainNodes)
            self.__nodalDomainDisp_Z = np.zeros(self.__nDomainNodes)
            self.__nodalDomainDispPrev_X = np.zeros(self.__nDomainNodes)
            self.__nodalDomainDispPrev_Y = np.zeros(self.__nDomainNodes)
            self.__nodalDomainDispPrev_Z = np.zeros(self.__nDomainNodes)
            self.__nodalDomainDispNm2_X = np.zeros(self.__nDomainNodes)
            self.__nodalDomainDispNm2_Y = np.zeros(self.__nDomainNodes)
            self.__nodalDomainDispNm2_Z = np.zeros(self.__nDomainNodes)
            self.__nodalDomainDispPrev_X, self.__nodalDomainDispPrev_Y, self.__nodalDomainDispPrev_Z = self.__vecToVecArray(GetDPGetNumber("nodalDisplacementPrev"))
            self.__nodalDomainDispNm2_X, self.__nodalDomainDispNm2_Y, self.__nodalDomainDispNm2_Z = self.__vecToVecArray(GetDPGetNumber("nodalDisplacementNm2"))
            self.__nodalDomainTemp = np.zeros(self.__nDomainNodes)
            self.__nodalDomainTempPrev = np.zeros(self.__nDomainNodes)
            self.__nodalDomainTempNm2 = np.zeros(self.__nDomainNodes)
            if self.regime == 'unsteady':
                self.__nodalDomainTemp = self.__vecToScalArray(GetDPGetNumber("nodalTemperatureNm0"))
                self.__nodalDomainTempPrev = self.__vecToScalArray(GetDPGetNumber("nodalTemperaturePrev"))
                self.__nodalDomainTempNm2 = self.__vecToScalArray(GetDPGetNumber("nodalTemperatureNm2"))
            self.nodalLoads_X = np.zeros(self.nPhysicalNodes)
            self.nodalLoads_Y = np.zeros(self.nPhysicalNodes)
            self.nodalLoads_Z = np.zeros(self.nPhysicalNodes) 
            self.__setCurrentState(True)
        else:
            os.system(self.pathToGetDP +" {} -setnumber Initialize 1 -setnumber OutputFiles 1 -solve {}".format(self.testname, self.resolution))
            with open("nodalPosition.txt") as f:
                strings = f.readlines()
            self.nNodes = len(strings)-1
            self.nHaloNode = 0
            self.nPhysicalNodes = len(strings)-1
            with open("nodalDisplacement.txt") as f:
                strings = f.readlines()
            self.__nDomainNodes = len(strings)-1
            SolidSolver.__init__(self)
            self.nodalInterfIndex = self.__readIndex("nodalPosition.txt")
            self.nodalDomainIndex = self.__readIndex("nodalDisplacement.txt")
            self.nodalInitialPos_X = np.zeros(self.nPhysicalNodes)
            self.nodalInitialPos_Y = np.zeros(self.nPhysicalNodes)
            self.nodalInitialPos_Z = np.zeros(self.nPhysicalNodes)
            self.nodalInitialPos_X, self.nodalInitialPos_Y, self.nodalInitialPos_Z = self.__readFileToVec("nodalPosition.txt", self.nPhysicalNodes)
            self.__nodalDomainDisp_X = np.zeros(self.__nDomainNodes)
            self.__nodalDomainDisp_Y = np.zeros(self.__nDomainNodes)
            self.__nodalDomainDisp_Z = np.zeros(self.__nDomainNodes)
            self.__nodalDomainDispPrev_X = np.zeros(self.__nDomainNodes)
            self.__nodalDomainDispPrev_Y = np.zeros(self.__nDomainNodes)
            self.__nodalDomainDispPrev_Z = np.zeros(self.__nDomainNodes)
            self.__nodalDomainDispNm2_X = np.zeros(self.__nDomainNodes)
            self.__nodalDomainDispNm2_Y = np.zeros(self.__nDomainNodes)
            self.__nodalDomainDispNm2_Z = np.zeros(self.__nDomainNodes)
            self.__nodalDomainDispPrev_X, self.__nodalDomainDispPrev_Y, self.__nodalDomainDispPrev_Z = self.__readFileToVec("nodalDisplacementPrev.txt", self.__nDomainNodes)
            self.__nodalDomainDispNm2_X, self.__nodalDomainDispNm2_Y, self.__nodalDomainDispNm2_Z = self.__readFileToVec("nodalDisplacementNm2.txt", self.__nDomainNodes)
            self.__nodalDomainTemp = np.zeros(self.__nDomainNodes)
            self.__nodalDomainTempPrev = np.zeros(self.__nDomainNodes)
            self.__nodalDomainTempNm2 = np.zeros(self.__nDomainNodes)
            if self.regime == 'unsteady':
                self.__nodalDomainTemp = self.__readFileToScal("nodalTemperatureNm0.txt", self.__nDomainNodes)
                self.__nodalDomainTempPrev = self.__readFileToScal("nodalTemperaturePrev.txt", self.__nDomainNodes)
                self.__nodalDomainTempNm2 = self.__readFileToScal("nodalTemperatureNm2.txt", self.__nDomainNodes)
            self.nodalLoads_X = np.zeros(self.nPhysicalNodes)
            self.nodalLoads_Y = np.zeros(self.nPhysicalNodes)
            self.nodalLoads_Z = np.zeros(self.nPhysicalNodes)
            self.__writeVecToFile("nodalHeatFlux.txt", np.zeros(self.nPhysicalNodes),  np.zeros(self.nPhysicalNodes), np.zeros(self.nPhysicalNodes), self.nodalInterfIndex)
            self.__setCurrentState(True)

        self.initRealTimeData()

    def __readIndex(self, fileName):


        dic = {}

        with open(fileName) as f:
            strings = f.readlines()

        iVertex = 0
        for iter in range(1, len(strings)):
            no = int(strings[iter].split()[0])
            dic[iVertex] = no
            iVertex += 1

        return dic

    def __extractIndex(self, vector, nDim):


        extractList = list(vector)
        indexList = extractList[1:(len(extractList)):nDim+1]
        dic = {}

        for iVertex in range(len(indexList)):
            no = int(indexList[iVertex])
            dic[iVertex] = no

        return dic

    def __readFileToVec(self, fileName, nNodes):


        vec_x = np.zeros(nNodes)
        vec_y = np.zeros(nNodes)
        vec_z = np.zeros(nNodes)

        with open(fileName) as f:
            strings = f.readlines()

        iVertex = 0
        for iter in range(1, len(strings)):
            no = int(strings[iter].split()[0])
            vec_x[iVertex] = float(strings[iter].split()[1])
            vec_y[iVertex] = float(strings[iter].split()[2])
            vec_z[iVertex] = float(strings[iter].split()[3])
            iVertex += 1

        return (vec_x, vec_y, vec_z)

    def __readFileToScal(self, fileName, nNodes):


        scal = np.zeros(nNodes)

        with open(fileName) as f:
            strings = f.readlines()

        iVertex = 0
        for iter in range(1, len(strings)):
            no = int(strings[iter].split()[0])
            scal[iVertex] = float(strings[iter].split()[1])
            iVertex += 1

        return scal

    def __vecToVecArray(self, vector):


        extractList = list(vector)
        nNodes = int(extractList[0])

        vec_x = np.array(extractList[2:(len(extractList)):4])
        vec_y = np.array(extractList[3:(len(extractList)):4])
        vec_z = np.array(extractList[4:(len(extractList)):4])

        return (vec_x, vec_y, vec_z)

    def __vecToScalArray(self, vector):


        extractList = list(vector)
        nNodes = int(extractList[0])

        vec = np.array(extractList[2:(len(extractList)):2])

        return vec

    def __writeVecToFile(self, fileName, vec_x, vec_y, vec_z, index):


        f = open(fileName, "w")
        nNodes = vec_x.size
        f.write(str(nNodes) + '\n')
        for iVertex in range(nNodes):
            node = index[iVertex]
            f.write(str(node) + " "+ str(vec_x[iVertex]) + " " + str(vec_y[iVertex]) + " " + str(vec_z[iVertex]) + "\n")
        f.close()

    def __writeScalToFile(self, fileName, scal, index):


        f = open(fileName, "w")
        nNodes = scal.size
        f.write(str(nNodes) + '\n')
        for iVertex in range(nNodes):
            node = index[iVertex]
            f.write(str(node) + " " + str(scal[iVertex]) + "\n" )
        f.close()

    def __vecArrayToVec(self, vec_x, vec_y, vec_z, index):


        nNodes = vec_x.size

        vec = list()
        vec.append(float(nNodes))
        for ii in range(nNodes):
            no = index[ii]
            vec.append(float(no))
            vec.append(vec_x[ii])
            vec.append(vec_y[ii])
            vec.append(vec_z[ii])

        return vec

    def __scalArrayToVec(self, scal, index):


        nNodes = scal.size

        vec = list()
        vec.append(float(nNodes))
        for ii in range(nNodes):
            no = index[ii]
            vec.append(float(no))       #really needs a float here ?
            vec.append(scal[ii])

        return vec

    def run(self, t1, t2):


        self.currentDT = t2-t1

        if self.pythonFlag:
            GetDPSetNumber("Initialize", 0)
            GetDPSetNumber("OutputFiles", 1)
            GetDPSetNumber("nodalDisplacementPrev", self.__vecArrayToVec(self.__nodalDomainDispPrev_X, self.__nodalDomainDispPrev_Y, self.__nodalDomainDispPrev_Z, self.nodalDomainIndex))
            GetDPSetNumber("nodalDisplacementNm2", self.__vecArrayToVec(self.__nodalDomainDispNm2_X, self.__nodalDomainDispNm2_Y, self.__nodalDomainDispNm2_Z, self.nodalDomainIndex))
            GetDPSetNumber("nodalTemperaturePrev", self.__scalArrayToVec(self.__nodalDomainTempPrev, self.nodalDomainIndex))
            GetDPSetNumber("nodalTemperatureNm2", self.__scalArrayToVec(self.__nodalDomainTempNm2, self.nodalDomainIndex))
            if self.regime == 'unsteady':
                GetDPSetNumber("T1", t1)
                GetDPSetNumber("T2", t2)
            GetDP(["getdp", self.testname, "-solve", self.resolution])
            self.__setCurrentState(False)
        else:
            self.__writeVecToFile("nodalDisplacementPrev.txt", self.__nodalDomainDispPrev_X, self.__nodalDomainDispPrev_Y, self.__nodalDomainDispPrev_Z, self.nodalDomainIndex)
            self.__writeVecToFile("nodalDisplacementNm2.txt", self.__nodalDomainDispNm2_X, self.__nodalDomainDispNm2_Y, self.__nodalDomainDispNm2_Z, self.nodalDomainIndex)
            self.__writeScalToFile("nodalTemperaturePrev.txt", self.__nodalDomainTempPrev, self.nodalDomainIndex)
            self.__writeScalToFile("nodalTemperatureNm2.txt", self.__nodalDomainTempNm2, self.nodalDomainIndex)
            if self.regime == 'unsteady':
                os.system(self.pathToGetDP +" {} -setnumber Initialize 0 -setnumber OutputFiles 1 -setnumber T1 {} -setnumber T2 {} -solve {}".format(self.testname, t1, t2, self.resolution))
            else:
                os.system(self.pathToGetDP +" {} -setnumber Initialize 0 -setnumber OutputFiles 1 -solve {}".format(self.testname, self.resolution))
            self.__setCurrentState(False)
        
        return True
            

    def __setCurrentState(self, initialize):


        if self.pythonFlag:
            nodalPos_X , nodalPos_Y, nodalPos_Z = self.__vecToVecArray(GetDPGetNumber("nodalPosition"))
            self.nodalVel_X, self.nodalVel_Y, self.nodalVel_Z = self.__vecToVecArray(GetDPGetNumber("nodalVelocity"))
            if self.regime == 'unsteady':
                self.__nodalDomainDisp_X, self.__nodalDomainDisp_Y, self.__nodalDomainDisp_Z = self.__vecToVecArray(GetDPGetNumber("nodalDisplacement"))
                self.__nodalDomainTemp = self.__vecToScalArray(GetDPGetNumber("nodalTemperatureNm0"))
            if not initialize:
                self.nodalHeatFlux_X , self.nodalHeatFlux_Y , self.nodalHeatFlux_Z = self.__vecToVecArray(GetDPGetNumber("nodalHeatFlux"))
                self.nodalTemperature = self.__vecToScalArray(GetDPGetNumber("nodalTemperature"))
        else:
            nodalPos_X , nodalPos_Y, nodalPos_Z = self.__readFileToVec("nodalPosition.txt", self.nPhysicalNodes)
            self.nodalVel_X, self.nodalVel_Y, self.nodalVel_Z = self.__readFileToVec("nodalVelocity.txt", self.nPhysicalNodes)
            if self.regime == 'unsteady':
                self.__nodalDomainDisp_X, self.__nodalDomainDisp_Y, self.__nodalDomainDisp_Z = self.__readFileToVec("nodalDisplacement.txt", self.__nDomainNodes)
                self.__nodalDomainTemp = self.__readFileToScal("nodalTemperatureNm0.txt", self.__nDomainNodes)
            self.nodalHeatFlux_X , self.nodalHeatFlux_Y , self.nodalHeatFlux_Z = self.__readFileToVec("nodalHeatFlux.txt", self.nPhysicalNodes)
            self.nodalTemperature = self.__readFileToScal("nodalTemperature.txt", self.nPhysicalNodes)

        self.nodalDisp_X = nodalPos_X - self.nodalInitialPos_X
        self.nodalDisp_Y = nodalPos_Y - self.nodalInitialPos_Y
        self.nodalDisp_Z = nodalPos_Z - self.nodalInitialPos_Z

    def getNodalInitialPositions(self):


        return (self.nodalInitialPos_X, self.nodalInitialPos_Y, self.nodalInitialPos_Z)

    def getNodalIndex(self, iVertex):


        return self.nodalInterfIndex[iVertex]

    def applyNodalForce(self, load_X, load_Y, load_Z, dt, haloNodesLoads):


        if self.pythonFlag:
            GetDPSetNumber("nodalForce", self.__vecArrayToVec(load_X, load_Y, load_Z, self.nodalInterfIndex))
        else:
            self.__writeVecToFile("nodalForce.txt", load_X, load_Y, load_Z, self.nodalInterfIndex)


    def applyNodalTemperatures(self, Temperature, dt, haloNodesTemperature):


        if self.pythonFlag:
            GetDPSetNumber("nodalTemperature", self.__scalArrayToVec(Temperature, self.nodalInterfIndex))
        else:
            self.__writeScalToFile("nodalTemperature.txt" , Temperature, self.nodalInterfIndex)

    def applyNodalNormalHeatFluxes(self, NormalHeatFlux, dt):


        if self.pythonFlag:
            GetDPSetNumber("nodalNormalHeatFlux", self.__scalArrayToVec(NormalHeatFlux, self.nodalInterfIndex))
        else:
            self.__writeScalToFile("nodalNormalHeatFlux.txt", NormalHeatFlux, self.nodalInterfIndex)

    def applyNodalHeatFluxes(self, HeatFlux_X, HeatFlux_Y, HeatFlux_Z, dt, haloNodesHeatFlux):


        if self.pythonFlag:
            GetDPSetNumber("nodalHeatFlux", self.__vecArrayToVec(HeatFlux_X, HeatFlux_Y, HeatFlux_Z, self.nodalInterfIndex))
        else:
            self.__writeVecToFile("nodalHeatFlux.txt", HeatFlux_X, HeatFlux_Y, HeatFlux_Z, self.nodalInterfIndex)

    def update(self):


        SolidSolver.update(self)

        self.__nodalDomainDispNm2_X = self.__nodalDomainDispPrev_X.copy()
        self.__nodalDomainDispNm2_Y = self.__nodalDomainDispPrev_Y.copy()
        self.__nodalDomainDispNm2_Z = self.__nodalDomainDispPrev_Z.copy()
        self.__nodalDomainDispPrev_X = self.__nodalDomainDisp_X.copy()
        self.__nodalDomainDispPrev_Y = self.__nodalDomainDisp_Y.copy()
        self.__nodalDomainDispPrev_Z = self.__nodalDomainDisp_Z.copy()

        self.__nodalDomainTempNm2 = self.__nodalDomainTempPrev.copy()
        self.__nodalDomainTempPrev = self.__nodalDomainTemp.copy()

    def initRealTimeData(self):


        self.extractNode = 7        #should be decided by the user outside of the function
        self.iVertexExtract = list(self.nodalInterfIndex.keys())[list(self.nodalInterfIndex.values()).index(self.extractNode)]
        solFile = open('extractPhysicalNode' + str(self.extractNode) + '.ascii', "w")
        solFile.write("Time\tnIter\tDx\tDy\tDz\tVx\tVy\tVz\tT\tQx\tQy\tQz\n")
        solFile.close()

    def saveRealTimeData(self, time, nFSIIter):


        solFile = open('extractPhysicalNode' + str(self.extractNode) + '.ascii', "a")
        
        Dx = self.nodalDisp_X[self.iVertexExtract]
        Dy = self.nodalDisp_Y[self.iVertexExtract]
        Dz = self.nodalDisp_Z[self.iVertexExtract]
        Vx = self.nodalVel_X[self.iVertexExtract]
        Vy = self.nodalVel_Y[self.iVertexExtract]
        Vz = self.nodalVel_Z[self.iVertexExtract]
        T = self.nodalTemperature[self.iVertexExtract]
        Qx = self.nodalHeatFlux_X[self.iVertexExtract]
        Qy = self.nodalHeatFlux_Y[self.iVertexExtract]
        Qz = self.nodalHeatFlux_Z[self.iVertexExtract]

        solFile.write(str(time) + '\t' + str(nFSIIter) + '\t' + str(Dx) + '\t' + str(Dy) + '\t' + str(Dz) + '\t' + str(Vx) + '\t' + str(Vy) + '\t' + str(Vz) + '\t' + str(T) + '\t' + str(Qx) + '\t' + str(Qy) + '\t' + str(Qz) + '\n')

        solFile.close()

    def exit(self):


        print("***************************** Exit GetDP *****************************")
