#!/usr/bin/env python
# -*- coding: latin-1; -*-

# \file GetDPSolver.py
#  \brief Python wrapping script used to sequentially run GetDP and exchange data with external fluid solver for coupled FSI computation.
#  \author C. GEUZAINE, D. THOMAS, University of Liege, Belgium. Department of Mechanical and Aerospace Engineering
#  \version BETA

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import math
import os
import numpy as np
from FSICoupler import SolidSolver

# ----------------------------------------------------------------------
#  GetDPSolver class
# ----------------------------------------------------------------------

class GetDPSolver(SolidSolver):
    """
    Des.
    """

    def __init__(self, testname, resolution, computationType, pythonFlag):
        """
        Des
        """

        self.testname = testname
        self.pythonFlag = pythonFlag
        self.computationType = computationType
        self.resolution = resolution
        self.forceN = {}
        self.dispN = {}
        self.dispNm1 = {}
        self.dispNm2 = {}
        self.stateN = {}
        self.stateNm1 = {}
        self.currentDT = 1.0
        self.pathToGetDP = "/home/dthomas/InstalledSoftware/GetDP/bin/getdp"

        if self.pythonFlag:
            import getdp
            GetDPSetNumber("Initialize", 1)
            GetDP(["getdp", self.testname, "-solve", self.resolution])
            p = GetDPGetNumber("nodalPosition");
            self.dispNm2 = self.__vecToDic(GetDPGetNumber("nodalDisplacementNm2"))
            self.dispNm1 = self.__vecToDic(GetDPGetNumber("nodalDisplacementNm1"))
            self.dispN = self.__vecToDic(GetDPGetNumber("nodalDisplacement"))  
            initPos =  self.__vecToDic(GetDPGetNumber("nodalPosition"))
            initVel =  self.__vecToDic(GetDPGetNumber("nodalVelocity"))
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
            self.__nodalDomainDispNm1_X = np.zeros(self.__nDomainNodes)
            self.__nodalDomainDispNm1_Y = np.zeros(self.__nDomainNodes)
            self.__nodalDomainDispNm1_Z = np.zeros(self.__nDomainNodes)
            self.__nodalDomainDispNm2_X = np.zeros(self.__nDomainNodes)
            self.__nodalDomainDispNm2_Y = np.zeros(self.__nDomainNodes)
            self.__nodalDomainDispNm2_Z = np.zeros(self.__nDomainNodes)
            self.__nodalDomainDispNm1_X, self.__nodalDomainDispNm1_Y, self.__nodalDomainDispNm1_Z = self.__readFileToVec("nodalDisplacementNm1.txt", self.__nDomainNodes)
            self.__nodalDomainDispNm2_X, self.__nodalDomainDispNm2_Y, self.__nodalDomainDispNm2_Z = self.__readFileToVec("nodalDisplacementNm2.txt", self.__nDomainNodes)
            self.__nodalDomainTemp = np.zeros(self.__nDomainNodes)
            self.__nodalDomainTempNm1 = np.zeros(self.__nDomainNodes)
            self.__nodalDomainTempNm2 = np.zeros(self.__nDomainNodes)
            self.__nodalDomainTemp = self.__readFileToScal("nodalTemperatureNm0.txt", self.__nDomainNodes)
            self.__nodalDomainTempNm1 = self.__readFileToScal("nodalTemperatureNm1.txt", self.__nDomainNodes)
            self.__nodalDomainTempNm2 = self.__readFileToScal("nodalTemperatureNm2.txt", self.__nDomainNodes)
            self.nodalLoads_X = np.zeros(self.nPhysicalNodes)
            self.nodalLoads_Y = np.zeros(self.nPhysicalNodes)
            self.nodalLoads_Z = np.zeros(self.nPhysicalNodes)
            self.__writeVecToFile("nodalHeatFlux.txt", np.zeros(self.nPhysicalNodes),  np.zeros(self.nPhysicalNodes), np.zeros(self.nPhysicalNodes), self.nodalInterfIndex)
            self.__setCurrentState()
            self.nodalVel_XNm1 = self.nodalVel_X.copy()
            self.nodalVel_YNm1 = self.nodalVel_Y.copy()
            self.nodalVel_ZNm1 = self.nodalVel_Z.copy()

        self.initRealTimeData()

    def __readIndex(self, fileName):
        """
        Des.
        """

        dic = {}

        with open(fileName) as f:
            strings = f.readlines()

        iVertex = 0
        for iter in range(1, len(strings)):
            no = int(strings[iter].split()[0])
            dic[iVertex] = no
            iVertex += 1

        return dic

    def __readFileToVec(self, fileName, nNodes):
        """
        Des.
        """

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
        """
        Des
        """

        scal = np.zeros(nNodes)

        with open(fileName) as f:
            strings = f.readlines()

        iVertex = 0
        for iter in range(1, len(strings)):
            no = int(strings[iter].split()[0])
            scal[iVertex] = float(strings[iter].split()[1])
            iVertex += 1

        return scal

    def __writeVecToFile(self, fileName, vec_x, vec_y, vec_z, index):
        """
        Des.
        """

        f = open(fileName, "w")
        nNodes = vec_x.size
        f.write(str(nNodes) + '\n')
        for iVertex in range(nNodes):
            node = index[iVertex]
            f.write(str(node) + " "+ str(vec_x[iVertex]) + " " + str(vec_y[iVertex]) + " " + str(vec_z[iVertex]) + "\n")
        f.close()

    def __writeScalToFile(self, fileName, scal, index):
        """
        Des.
        """

        f = open(fileName, "w")
        nNodes = scal.size
        f.write(str(nNodes) + '\n')
        for iVertex in range(nNodes):
            node = index[iVertex]
            f.write(str(node) + " " + str(scal[iVertex]) + "\n" )
        f.close()

    def run(self, t1, t2):
        """
        Des.
        """

        self.currentDT = t2-t1

        if self.pythonFlag:
            hi=0
        else:
            self.__writeVecToFile("nodalDisplacementNm1.txt", self.__nodalDomainDispNm1_X, self.__nodalDomainDispNm1_Y, self.__nodalDomainDispNm1_Z, self.nodalDomainIndex)
            self.__writeVecToFile("nodalDisplacementNm2.txt", self.__nodalDomainDispNm2_X, self.__nodalDomainDispNm2_Y, self.__nodalDomainDispNm2_Z, self.nodalDomainIndex)
            self.__writeScalToFile("nodalTemperatureNm1.txt", self.__nodalDomainTempNm1, self.nodalDomainIndex)
            self.__writeScalToFile("nodalTemperatureNm2.txt", self.__nodalDomainTempNm2, self.nodalDomainIndex)
            if self.computationType == 'unsteady':
                os.system(self.pathToGetDP +" {} -setnumber Initialize 0 -setnumber OutputFiles 1 -setnumber T1 {} -setnumber T2 {} -solve {}".format(self.testname, t1, t2, self.resolution))
            else:
                os.system(self.pathToGetDP +" {} -setnumber Initialize 0 -setnumber OutputFiles 1 -solve {}".format(self.testname, self.resolution))
            self.__setCurrentState()
            

    def __setCurrentState(self):
        """
        Des.
        """

        nodalPos_X , nodalPos_Y, nodalPos_Z = self.__readFileToVec("nodalPosition.txt", self.nPhysicalNodes)

        self.nodalDisp_X = nodalPos_X - self.nodalInitialPos_X
        self.nodalDisp_Y = nodalPos_Y - self.nodalInitialPos_Y
        self.nodalDisp_Z = nodalPos_Z - self.nodalInitialPos_Z

        self.nodalVel_X, self.nodalVel_Y, self.nodalVel_Z = self.__readFileToVec("nodalVelocity.txt", self.nPhysicalNodes)

        self.__nodalDomainDisp_X, self.__nodalDomainDisp_Y, self.__nodalDomainDisp_Z = self.__readFileToVec("nodalDisplacement.txt", self.__nDomainNodes)

        self.__nodalDomainTemp = self.__readFileToScal("nodalTemperatureNm0.txt", self.__nDomainNodes)

        self.nodalHeatFlux_X , self.nodalHeatFlux_Y , self.nodalHeatFlux_Z = self.__readFileToVec("nodalHeatFlux.txt", self.nPhysicalNodes)

        self.nodalTemperature = self.__readFileToScal("nodalTemperature.txt", self.nPhysicalNodes)

    def getNodalInitialPositions(self):
        """
        Des.
        """

        return (self.nodalInitialPos_X, self.nodalInitialPos_Y, self.nodalInitialPos_Z)

    def getNodalIndex(self, iVertex):
        """
        des.
        """

        return self.nodalInterfIndex[iVertex]

    def applyNodalLoads(self, load_X, load_Y, load_Z, time):
        """
        Des.
        """

        self.__writeVecToFile("nodalForce.txt", load_X, load_Y, load_Z, self.nodalInterfIndex)


    def applyNodalTemperatures(self, Temperature, val_time):
        """
        Des.
        """

        self.__writeScalToFile("nodalTemperature.txt" , Temperature, self.nodalInterfIndex)

    def applyNodalNormalHeatFluxes(self, NormalHeatFlux, val_time):
        """
        Des.
        """

        self.__writeScalToFile("nodalNormalHeatFlux.txt", NormalHeatFlux, self.nodalInterfIndex)

    def applyNodalHeatFluxes(self, HeatFlux_X, HeatFlux_Y, HeatFlux_Z, time):
        """
        Des.
        """

        self.__writeVecToFile("nodalHeatFlux.txt", HeatFlux_X, HeatFlux_Y, HeatFlux_Z)

    def update(self):
        """
        Des.
        """

        SolidSolver.update(self)

        self.__nodalDomainDispNm2_X = self.__nodalDomainDispNm1_X.copy()
        self.__nodalDomainDispNm2_Y = self.__nodalDomainDispNm1_Y.copy()
        self.__nodalDomainDispNm2_Z = self.__nodalDomainDispNm1_Z.copy()
        self.__nodalDomainDispNm1_X = self.__nodalDomainDisp_X.copy()
        self.__nodalDomainDispNm1_Y = self.__nodalDomainDisp_Y.copy()
        self.__nodalDomainDispNm1_Z = self.__nodalDomainDisp_Z.copy()

        self.__nodalDomainTempNm2 = self.__nodalDomainTempNm1.copy()
        self.__nodalDomainTempNm1 = self.__nodalDomainTemp.copy()

    def initRealTimeData(self):
        """
        des.
        """

        self.extractNode = 7        #should be decided by the user outside of the function
        self.iVertexExtract = self.nodalInterfIndex.keys()[self.nodalInterfIndex.values().index(self.extractNode)]
        solFile = open('extractPhysicalNode' + str(self.extractNode) + '.ascii', "w")
        solFile.write("Time\tnIter\tDx\tDy\tDz\tVx\tVy\tVz\tT\tQx\tQy\tQz\n")
        solFile.close()

    def saveRealTimeData(self, time, nFSIIter):
        """
        des.
        """

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
        """
        Des.
        """

        print("***************************** Exit GetDP *****************************")
