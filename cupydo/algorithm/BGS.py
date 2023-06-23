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

algorithm.py
Defines the coupling algorithms of CUPyDO.
Authors : David THOMAS, Marco Lucio CERQUAGLIA, Romain BOMAN

'''

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

from math import *
import numpy as np
import sys
import traceback
import sys

from ..utilities import *
from ..interfaceData import FlexInterfaceData
from .algorithm import Algorithm

np.set_printoptions(threshold=sys.maxsize)

# ----------------------------------------------------------------------
#    Static BGS Algorithm class
# ----------------------------------------------------------------------

class AlgorithmBGSStaticRelax(Algorithm):

    def __init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, Criterion, nbFSIIterMax, deltaT, totTime, dtSave, omegaBoundList=[1.0,1.0], mpiComm=None):

        Algorithm.__init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, deltaT, totTime, dtSave, mpiComm)

        self.criterion = Criterion

        if type(omegaBoundList) == list:        #A list specified by the user or the default one
            self.omegaBoundMecha = omegaBoundList[0]
            self.omegaBoundThermal = omegaBoundList[1]
        else:                                   #The user just put one value
            self.omegaBoundMecha = omegaBoundList
            self.omegaBoundThermal = 1.0
        self.omegaMinMecha = 1e-12
        self.omegaMinThermal = 1e-12
        self.omegaMecha = self.omegaBoundMecha
        self.omegaThermal = self.omegaBoundThermal

        self.writeInFSIloop = False

        self.FSIIter = 0
        self.errValue = 0.0
        self.FSIConv = False
        self.totNbOfFSIIt = 0
        self.nbFSIIterMax = nbFSIIterMax

        self.predictor = True
        self.predictorOrder = 2
        self.alpha_0 = 1.0
        self.alpha_1 = 0.5

        self.solidInterfaceVelocity = None
        self.solidInterfaceVelocitynM1 = None
        self.solidInterfaceResidual = None
        self.solidHeatFluxResidual = None
        self.solidTemperatureResidual = None

    def initInterfaceData(self):

        ns = self.interfaceInterpolator.getNs()
        d = self.interfaceInterpolator.getd()

        # --- Initialize data for prediction (mechanical only) --- #
        if self.predictor and self.manager.mechanical:
            self.solidInterfaceVelocity = FlexInterfaceData(ns+d, 3, self.mpiComm)
            self.solidInterfaceVelocitynM1 = FlexInterfaceData(ns+d, 3, self.mpiComm)

        # --- Initialize coupling residuals --- #
        if self.manager.mechanical:
            self.solidInterfaceResidual = FlexInterfaceData(ns+d, 3, self.mpiComm)
        if self.manager.thermal:
            self.solidHeatFluxResidual = FlexInterfaceData(ns+d, 3, self.mpiComm)
            self.solidTemperatureResidual = FlexInterfaceData(ns+d, 1, self.mpiComm)

    def run(self):

        # --- Initialize the algorithm --- #
        mpiPrint("Begin FSI Computation",self.mpiComm,titlePrint)
        self.initInterfaceData()
        self.iniRealTimeData()
        self.globalTimer.start()
        self.setFSIInitialConditions()

        try:
            if self.manager.computationType == 'unsteady':
                self.__unsteadyRun()
            else:
                self.step.dt = self.totTime
                self.writeInFSIloop = True

                # --- Internal FSI loop --- #
                self.verified = self.fsiCoupling()
                self.totNbOfFSIIt = self.FSIIter
                if not self.verified: raise Exception('The steady FSI coupling did not converge')

                self.FluidSolver.save(self.step.timeIter)
                if self.myid in self.manager.getSolidSolverProcessors():
                    self.SolidSolver.save()
    
        except:
            mpiPrint('\nError when executing BGS : unsteadyRun\n', self.mpiComm)
            traceback.print_exc()
        finally:
            self.globalTimer.stop()
            self.globalTimer.cumul()
            
            mpiBarrier(self.mpiComm)
            mpiPrint("End FSI Computation",self.mpiComm,titlePrint)
            self.printExitInfo()
            
            # --- Exit the solid solver --- #
            if self.myid in self.manager.getSolidSolverProcessors():
                self.SolidSolver.exit()

            # --- Exit the fluid solver --- #
            self.FluidSolver.exit()
    
            # --- Exit computation --- #
            mpiBarrier(self.mpiComm)

    def __unsteadyRun(self):

        mpiPrint('Begin time integration\n', self.mpiComm)

        # --- External temporal loop --- #
        while self.step.time < self.totTime:
            
            mpiPrint("\n>>>> Time iteration {} <<<<".format(self.step.timeIter), self.mpiComm)

            # --- Displacement predictor for the next time step --- #
            mpiPrint('\nSolid displacement prediction for next time step', self.mpiComm)
            self.solidDisplacementPredictor()

            # --- Preprocess the temporal iteration --- #
            self.FluidSolver.preprocessTimeIter(self.step.timeIter)
            if self.myid in self.manager.getSolidSolverProcessors():
                self.SolidSolver.preprocessTimeIter(self.step.timeIter)

            # --- Internal FSI loop --- #
            self.verified = self.fsiCoupling()
            self.totNbOfFSIIt += self.FSIIter
            mpiBarrier(self.mpiComm)

            # --- Update TimeStep class and restart if FSI failed --- #
            if not self.verified:
                self.step.updateTime(self.verified)
                self.resetInternalVars()
                self.writeRealTimeData()
                continue

            # --- Update the fluid and solid solver for the next time step --- #
            if self.myid in self.manager.getSolidSolverProcessors():
                self.SolidSolver.update()
            self.FluidSolver.update(self.step.dt)

            # --- Perform some remeshing if necessary
            if self.myid in self.manager.getSolidSolverProcessors():
                self.solidRemeshingTimer.start()
                self.SolidSolver.remeshing()
                self.solidRemeshingTimer.stop()
                self.solidRemeshingTimer.cumul()
            
            self.fluidRemeshingTimer.start()
            self.FluidSolver.remeshing()
            self.fluidRemeshingTimer.stop()
            self.fluidRemeshingTimer.cumul()

            # --- Update TimeStep class, export the results and write FSI history --- #
            self.step.updateTime(self.verified)
            self.step.updateSave(self)
            self.writeRealTimeData()

    def iniRealTimeData(self):

        if self.myid == 0:
            self.FluidSolver.initRealTimeData()
        if self.myid in self.manager.getSolidSolverProcessors():
            self.SolidSolver.initRealTimeData()
        histFile = open('FSIhistory.ascii', "w")
        histFile.write("{0:>12s}   {1:>12s}   {2:>12s}   {3:>12s}   {4:>12s}   {5:>12s}   {6:>12s}\n".format("TimeIter", "Time", "FSIError", "CHTError", "FSINbIter", "omegaMecha", "omegaThermal"))
        histFile.close()

    def writeRealTimeData(self):

        if self.myid == 0:
            self.FluidSolver.saveRealTimeData(self.step.time, self.FSIIter)
            self.SolidSolver.saveRealTimeData(self.step.time, self.FSIIter)
            histFile = open('FSIhistory.ascii', "a")
            histFile.write('{0:12d}   {1:.6e}   {2:.6e}   {3:.6e}   {4:12d}   {5:.6e}   {6:.6e}\n'.format(self.step.timeIter, self.step.time, self.errValue, self.errValue_CHT, self.FSIIter, self.omegaMecha, self.omegaThermal))
            histFile.close()

    def getMeanNbOfFSIIt(self):

        if self.manager.computationType == 'unsteady':
            return float(self.totNbOfFSIIt)/(self.step.timeIter+1)
        else:
            return self.FSIIter

    def printExitInfo(self):

        mpiPrint('[cpu FSI total]: ' + str(self.globalTimer.cumulTime) + ' s', self.mpiComm)
        mpiPrint('[cpu FSI fluid mesh mapping]: ' + str(self.interfaceInterpolator.mappingTimer.cumulTime) + ' s', self.mpiComm)
        mpiPrint('[cpu FSI fluid mesh deformation]: ' + str(self.meshDefTimer.cumulTime) + ' s', self.mpiComm)
        mpiPrint('[cpu FSI communications]: ' + str(self.communicationTimer.cumulTime) + ' s', self.mpiComm)
        mpiPrint('[cpu FSI fluid solver]: ' + str(self.fluidSolverTimer.cumulTime) + ' s', self.mpiComm)
        mpiPrint('[cpu FSI solid solver]: ' + str(self.solidSolverTimer.cumulTime) + ' s', self.mpiComm)
        mpiPrint('[cpu FSI fluid remeshing]: ' + str(self.fluidRemeshingTimer.cumulTime) + ' s', self.mpiComm)
        mpiPrint('[cpu FSI solid remeshing]: ' + str(self.solidRemeshingTimer.cumulTime) + ' s', self.mpiComm)
        mpiPrint('[Time steps FSI]: ' + str(self.step.timeIter), self.mpiComm)
        mpiPrint('[Successful Run FSI]: ' + str(self.step.time >= (self.totTime - 2*self.step.dt)), self.mpiComm) # NB: self.totTime - 2*self.step.dt is the extreme case that can be encountered due to rounding effects!
        mpiPrint('[Mean n. of FSI Iterations]: ' + str(self.getMeanNbOfFSIIt()), self.mpiComm)

        if self.myid == 0 :
            self.FluidSolver.printRealTimeData(self.step.time, self.FSIIter)
            self.SolidSolver.printRealTimeData(self.step.time, self.FSIIter)

        mpiPrint('RES-FSI-FSIhistory: ' + str(self.step.timeIter) + '\t' + str(self.step.time) + '\t' + str(self.errValue) + '\t' + str(self.FSIIter) + '\n', self.mpiComm)

    def fsiCoupling(self):
        """
        Block Gauss Seidel (BGS) method for strong coupling FSI
        """

        mpiPrint("Enter BGS Strong Coupling FSI",self.mpiComm,titlePrint)

        verif = False
        self.FSIIter = 0
        self.errValue = 1e12
        self.errValue_CHT = 1e6

        while (self.FSIIter < self.nbFSIIterMax):
            mpiPrint("\n>>>> FSI iteration {} <<<<\n".format(self.FSIIter), self.mpiComm)

            if self.manager.mechanical:
                # --- Solid to fluid mechanical transfer --- #
                self.solidToFluidMechaTransfer()
                # --- Fluid mesh morphing --- #
                mpiPrint('\nPerforming mesh deformation...\n', self.mpiComm)
                self.meshDefTimer.start()
                self.FluidSolver.meshUpdate(self.step.timeIter)
                self.meshDefTimer.stop()
                self.meshDefTimer.cumul()
            if self.manager.thermal and self.solidHasRun:
                # --- Solid to fluid thermal transfer --- #
                self.solidToFluidThermalTransfer()
            self.FluidSolver.boundaryConditionsUpdate()

            # --- Fluid solver call for FSI subiteration --- #
            mpiPrint('\nLaunching fluid solver...', self.mpiComm)
            self.fluidSolverTimer.start()
            verif = self.FluidSolver.run(*self.step.timeFrame())
            self.fluidSolverTimer.stop()
            self.fluidSolverTimer.cumul()
            mpiBarrier(self.mpiComm)

            # --- Check if the fluid solver succeeded --- #
            if not verif: return False

            if self.manager.mechanical:
                # --- Fluid to solid mechanical transfer --- #
                mpiPrint('\nProcessing interface fluid loads...\n', self.mpiComm)
                self.fluidToSolidMechaTransfer()
            if self.manager.thermal:
                # --- Fluid to solid thermal transfer --- #
                mpiPrint('\nProcessing interface thermal quantities...\n', self.mpiComm)
                self.fluidToSolidThermalTransfer()
            mpiBarrier(self.mpiComm)

            # --- Solid solver call for FSI subiteration --- #
            mpiPrint('\nLaunching solid solver...\n', self.mpiComm)
            if self.myid in self.manager.getSolidSolverProcessors():
                self.solidSolverTimer.start()
                verif = self.SolidSolver.run(*self.step.timeFrame())
                self.solidSolverTimer.stop()
                self.solidSolverTimer.cumul()

            # --- Check if the solid solver succeeded --- #
            try: solidProc = int(self.manager.getSolidSolverProcessors())
            except: raise Exception('Only one solid solver process is supported yet')
            verif = mpiScatter(verif, self.mpiComm, solidProc)
            self.solidHasRun = True
            if not verif: return False

            if self.manager.mechanical:
                # --- Compute the mechanical residual --- #
                res = self.computeSolidInterfaceResidual()
                self.errValue = self.criterion.update(res)
                mpiPrint('\nFSI error value : {}\n'.format(self.errValue), self.mpiComm)
            else:
                self.errValue = 0.0
            if self.manager.thermal:
                # --- Compute the thermal residual --- #
                res_CHT = self.computeSolidInterfaceResidual_CHT()
                self.errValue_CHT = self.criterion.updateThermal(res_CHT)
                mpiPrint('\nCHT error value : {}\n'.format(self.errValue_CHT), self.mpiComm)
            else:
                self.errValue_CHT = 0.0

            if self.manager.mechanical:
                # --- Relaxe the solid position --- #
                mpiPrint('\nProcessing interface displacements...\n', self.mpiComm)
                self.relaxSolidPosition()

            if self.manager.thermal:
                # --- Relaxe thermal data --- #
                self.relaxCHT()

            # --- Update the solvers for the next BGS steady iteration --- #
            if self.manager.computationType == 'steady':

                if self.myid in self.manager.getSolidSolverProcessors():
                    self.SolidSolver.steadyUpdate()
                self.FluidSolver.steadyUpdate()

            # --- Update the FSI iteration and history --- #
            if self.writeInFSIloop == True:
                self.writeRealTimeData()
            self.FSIIter += 1

            # --- Monitor the coupling convergence --- #
            if self.criterion.isVerified(self.errValue, self.errValue_CHT):
                mpiPrint("BGS is Converged",self.mpiComm,titlePrint)
                return True
        
        return False

    def computeSolidInterfaceResidual(self):

        ns = self.interfaceInterpolator.getNs()
        d = self.interfaceInterpolator.getd()

        # --- Get the predicted (computed) solid interface displacement from the solid solver --- #
        predictedDisplacement = FlexInterfaceData(ns+d, 3, self.mpiComm)

        if self.myid in self.manager.getSolidInterfaceProcessors():
            localSolidInterfaceDisp_X, localSolidInterfaceDisp_Y, localSolidInterfaceDisp_Z = self.SolidSolver.getNodalDisplacements()
            for iVertex in range(self.manager.getNumberOfLocalSolidInterfaceNodes()):
                iGlobalVertex = self.manager.getGlobalIndex('solid', self.myid, iVertex)
                predictedDisplacement[iGlobalVertex] = [localSolidInterfaceDisp_X[iVertex], localSolidInterfaceDisp_Y[iVertex], localSolidInterfaceDisp_Z[iVertex]]

        predictedDisplacement.assemble()

        # --- Calculate the residual (vector and norm) --- #
        mpiPrint("\nCompute FSI residual based on solid interface displacement.", self.mpiComm)
        self.solidInterfaceResidual.set(predictedDisplacement - self.interfaceInterpolator.solidInterfaceDisplacement)

        return self.solidInterfaceResidual

    def computeSolidInterfaceResidual_CHT(self):

        ns = self.interfaceInterpolator.getNs()
        d = self.interfaceInterpolator.getd()

        predictedHF = FlexInterfaceData(ns+d, 3, self.mpiComm)
        predictedTemp = FlexInterfaceData(ns+d, 1, self.mpiComm)

        if self.myid in self.manager.getSolidInterfaceProcessors():
            localSolidInterfaceHeatFlux_X, localSolidInterfaceHeatFlux_Y, localSolidInterfaceHeatFlux_Z = self.SolidSolver.getNodalHeatFluxes()
            localSolidInterfaceTemperature = self.SolidSolver.getNodalTemperatures()
            for iVertex in range(self.manager.getNumberOfLocalSolidInterfaceNodes()):
                iGlobalVertex = self.manager.getGlobalIndex('solid', self.myid, iVertex)
                predictedHF[iGlobalVertex] = [localSolidInterfaceHeatFlux_X[iVertex], localSolidInterfaceHeatFlux_Y[iVertex], localSolidInterfaceHeatFlux_Z[iVertex]]
                predictedTemp[iGlobalVertex] = [localSolidInterfaceTemperature[iVertex]]

        predictedHF.assemble()
        predictedTemp.assemble()

        if self.interfaceInterpolator.chtTransferMethod == 'hFFB' or self.interfaceInterpolator.chtTransferMethod == 'TFFB':
            mpiPrint("\nCompute CHT residual based on solid interface heat flux.", self.mpiComm)
            self.solidHeatFluxResidual.set(predictedHF - self.interfaceInterpolator.solidInterfaceHeatFlux)
            return self.solidHeatFluxResidual
        elif self.interfaceInterpolator.chtTransferMethod == 'hFTB' or self.interfaceInterpolator.chtTransferMethod == 'FFTB':
            mpiPrint("\nCompute CHT residual based on solid interface temperature.", self.mpiComm)
            self.solidTemperatureResidual.set(predictedTemp - self.interfaceInterpolator.solidInterfaceTemperature)
            return self.solidTemperatureResidual
        else:
            return None

    def solidDisplacementPredictor(self):

        if not self.predictor:
            if not self.verified:

                # --- Bring back the interface position of the previous time step --- #
                self.interfaceInterpolator.restoreSolidInterfaceDisplacement()
            
            return

        if self.verified:
            
            # --- Save the current interface position before predicting a new one --- #
            self.interfaceInterpolator.saveSolidInterfaceDisplacement()

            # --- Get the velocity (current and previous time step) of the solid interface from the solid solver --- #
            if self.myid in self.manager.getSolidInterfaceProcessors():
                localSolidInterfaceVel_X, localSolidInterfaceVel_Y, localSolidInterfaceVel_Z = self.SolidSolver.getNodalVelocity()
                localSolidInterfaceVelNm1_X, localSolidInterfaceVelNm1_Y, localSolidInterfaceVelNm1_Z = self.SolidSolver.getNodalVelocityNm1()
                for iVertex in range(self.manager.getNumberOfLocalSolidInterfaceNodes()):
                    iGlobalVertex = self.manager.getGlobalIndex('solid', self.myid, iVertex)
                    self.solidInterfaceVelocity[iGlobalVertex] = [localSolidInterfaceVel_X[iVertex], localSolidInterfaceVel_Y[iVertex], localSolidInterfaceVel_Z[iVertex]]
                    self.solidInterfaceVelocitynM1[iGlobalVertex] = [localSolidInterfaceVelNm1_X[iVertex], localSolidInterfaceVelNm1_Y[iVertex], localSolidInterfaceVelNm1_Z[iVertex]]

            self.solidInterfaceVelocity.assemble()
            self.solidInterfaceVelocitynM1.assemble()

        else:

            # --- Bring back the interface position of the previous time step --- #
            self.interfaceInterpolator.restoreSolidInterfaceDisplacement()

        # --- Predict the solid position for the next time step --- #
        if self.predictorOrder == 1:

            mpiPrint("First order predictor.", self.mpiComm)
            self.interfaceInterpolator.solidInterfaceDisplacement += (self.alpha_0*self.step.dt*self.solidInterfaceVelocity)

        elif self.predictorOrder == 2:

            mpiPrint("Second order predictor.", self.mpiComm)
            self.interfaceInterpolator.solidInterfaceDisplacement += (self.alpha_0*self.step.dt*self.solidInterfaceVelocity + self.alpha_1*self.step.dt*(self.solidInterfaceVelocity-self.solidInterfaceVelocitynM1))

        else:
            raise Exception('Only first or second order prdictors are available')

    def setOmegaMecha(self):

        self.omegaMecha = self.omegaBoundMecha
        mpiPrint('Static under-relaxation summary, mechanical : {}'.format(self.omegaMecha), self.mpiComm)


    def setOmegaThermal(self):

        self.omegaThermal = self.omegaBoundThermal
        mpiPrint('Static under-relaxation summary, thermal : {}'.format(self.omegaThermal), self.mpiComm)

    def relaxSolidPosition(self):

        # --- Set the relaxation parameter --- #
        self.setOmegaMecha()

        # --- Relax the solid interface position --- #
        self.interfaceInterpolator.solidInterfaceDisplacement += (self.omegaMecha*self.solidInterfaceResidual)

    def relaxCHT(self):

        self.setOmegaThermal()

        if self.interfaceInterpolator.chtTransferMethod == 'hFFB' or self.interfaceInterpolator.chtTransferMethod == 'TFFB':
            self.interfaceInterpolator.solidInterfaceHeatFlux += (self.omegaThermal*self.solidHeatFluxResidual)
        elif self.interfaceInterpolator.chtTransferMethod == 'hFTB' or self.interfaceInterpolator.chtTransferMethod == 'FFTB':
            self.interfaceInterpolator.solidInterfaceTemperature += (self.omegaThermal*self.solidTemperatureResidual)

# ----------------------------------------------------------------------
#    Aitken BGS Algorithm class
# ----------------------------------------------------------------------

class AlgorithmBGSAitkenRelax(AlgorithmBGSStaticRelax):

    def __init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, Criterion, nbFSIIterMax, deltaT, totTime, dtSave, omegaBoundList=[1.0, 1.0], mpiComm=None):

        AlgorithmBGSStaticRelax.__init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, Criterion, nbFSIIterMax, deltaT, totTime, dtSave, omegaBoundList, mpiComm)

        self.solidInterfaceResidualkM1 = None
        self.solidHeatFluxResidualkM1 = None
        self.solidTemperatureResidualkM1 = None
        self.aitkenCritMecha = 'max'
        self.aitkenCritThermal = 'max'

    def initInterfaceData(self):

        AlgorithmBGSStaticRelax.initInterfaceData(self)
        ns = self.interfaceInterpolator.getNs()
        d = self.interfaceInterpolator.getd()

        self.solidInterfaceResidualkM1 = FlexInterfaceData(ns+d, 3, self.mpiComm)
        self.solidHeatFluxResidualkM1 = FlexInterfaceData(ns+d, 3, self.mpiComm)
        self.solidTemperatureResidualkM1 = FlexInterfaceData(ns+d, 1, self.mpiComm)

    def resetInternalVars(self):

        self.omegaMecha = self.omegaBoundMecha
        self.omegaThermal = self.omegaBoundThermal

    def setOmegaMecha(self):

        if self.FSIIter != 0:
            # --- Compute the dynamic Aitken coefficient --- #
            deltaInterfaceResidual = self.solidInterfaceResidual - self.solidInterfaceResidualkM1

            prodScalRes_X, prodScalRes_Y, prodScalRes_Z = deltaInterfaceResidual.dot(self.solidInterfaceResidualkM1)
            prodScalRes = prodScalRes_X + prodScalRes_Y + prodScalRes_Z

            deltaInterfaceResidual_NormX, deltaInterfaceResidual_NormY, deltaInterfaceResidual_NormZ = deltaInterfaceResidual.norm()
            deltaResNormSquare = deltaInterfaceResidual_NormX**2 + deltaInterfaceResidual_NormY**2 + deltaInterfaceResidual_NormZ**2

            if deltaResNormSquare != 0.:
                self.omegaMecha *= -prodScalRes/deltaResNormSquare
            else:
                self.omegaMecha = self.omegaMinMecha

        else:
            # --- Initiate omega with min/max bounding --- #
            if self.aitkenCritMecha == 'max':
                self.omegaMecha = max(self.omegaBoundMecha, self.omegaMecha)
            else:
                self.omega = min(self.omegaBoundMecha, self.omegaMecha)

        self.omegaMecha = min(self.omegaMecha, 1.0)
        self.omegaMecha = max(self.omegaMecha, self.omegaMinMecha)

        mpiPrint('Aitken under-relaxation summary, mechanical : {}'.format(self.omegaMecha), self.mpiComm)

        # --- Update the value of the residual for the next FSI iteration --- #
        self.solidInterfaceResidual.copy(self.solidInterfaceResidualkM1)

    def setOmegaThermal(self):

        if self.FSIIter != 0:
            # --- Compute the dynamic Aitken coefficient --- #
            if self.interfaceInterpolator.chtTransferMethod == 'hFFB' or self.interfaceInterpolator.chtTransferMethod == 'TFFB':
                deltaHeatFluxResidual = self.solidHeatFluxResidual - self.solidHeatFluxResidualkM1

                prodScalRes_X, prodScalRes_Y, prodScalRes_Z = deltaHeatFluxResidual.dot(self.solidHeatFluxResidualkM1)
                prodScalRes = prodScalRes_X + prodScalRes_Y + prodScalRes_Z

                deltaHeatFluxResidual_NormX, deltaHeatFluxResidual_NormY, deltaHeatFluxResidual_NormZ = deltaHeatFluxResidual.norm()
                deltaResNormSquare = deltaHeatFluxResidual_NormX**2 + deltaHeatFluxResidual_NormY**2 + deltaHeatFluxResidual_NormZ**2
            elif self.interfaceInterpolator.chtTransferMethod == 'hFTB' or self.interfaceInterpolator.chtTransferMethod == 'FFTB':
                deltaTemperatureResidual = self.solidTemperatureResidual - self.solidTemperatureResidualkM1

                tempDot = deltaTemperatureResidual.dot(self.solidTemperatureResidualkM1)
                prodScalRes = tempDot[0]

                tempNorm = deltaTemperatureResidual.norm()
                deltaTemperatureResidual_Norm = tempNorm[0]
                deltaResNormSquare = deltaTemperatureResidual_Norm**2

            if deltaResNormSquare != 0.:
                self.omegaThermal *= -prodScalRes/deltaResNormSquare
            else:
                self.omegaThermal = self.omegaMinThermal

        else:
            # --- Initiate omega with min/max bounding --- #
            if self.aitkenCritThermal == 'max':
                self.omegaThermal = max(self.omegaBoundThermal, self.omegaThermal)
            else:
                self.omegaThermal = min(self.omegaBoundThermal, self.omegaThermal)

        self.omegaThermal = min(self.omegaThermal, 1.0)
        self.omegaThermal = max(self.omegaThermal, self.omegaMinThermal)

        mpiPrint('Aitken under-relaxation summary, thermal : {}'.format(self.omegaThermal), self.mpiComm)

        # --- Update the value of the residual for the next FSI iteration --- #
        self.solidHeatFluxResidual.copy(self.solidInterfaceResidualkM1)
        self.solidTemperatureResidual.copy(self.solidTemperatureResidualkM1)

class AlgorithmBGSStaticRelaxAdjoint(AlgorithmBGSStaticRelax):
    def run(self):

        # --- Initialize the algorithm --- #
        mpiPrint("Begin FSI Computation",self.mpiComm,titlePrint)
        self.initInterfaceData()
        self.iniRealTimeData()
        self.globalTimer.start()
        self.setFSIInitialConditions()
        
        try:
            self.step.timeIter = 1
            self.step.dt = self.totTime
            self.writeInFSIloop = True

            # --- Internal FSI loop --- #
            self.verified = self.fsiCoupling()
            self.totNbOfFSIIt = self.FSIIter
            if not self.verified: raise Exception('The adjoint FSI coupling did not converge')

            self.FluidSolver.save(self.step.timeIter)
            if self.myid in self.manager.getSolidSolverProcessors():
                self.SolidSolver.save()
    
        except:
            mpiPrint('\nError when executing BGS Adjoint : fsiCoupling\n', self.mpiComm)
            traceback.print_exc()
        finally:
            self.globalTimer.stop()
            self.globalTimer.cumul()
            
            mpiBarrier(self.mpiComm)
            mpiPrint("End FSI Adjoint Computation",self.mpiComm,titlePrint)
            self.printExitInfo()
            
            # --- Exit the solid solver --- #
            if self.myid in self.manager.getSolidSolverProcessors():
                self.SolidSolver.exit()

            # --- Exit the fluid solver --- #
            self.FluidSolver.exit()
    
            # --- Exit computation --- #
            mpiBarrier(self.mpiComm)

    def fsiCoupling(self):
        """
        Block Gauss Seidel (BGS) method for strong coupling FSI
        """

        mpiPrint("Enter BGS Strong Coupling FSI",self.mpiComm,titlePrint)

        self.FSIIter = 0
        self.FSIConv = False
        self.errValue = 1e12
        self.errValue_CHT = 1e6

        while ((self.FSIIter < self.nbFSIIterMax) and (not self.criterion.isVerified(self.errValue, self.errValue_CHT))):
            mpiPrint("\n>>>> FSI Adjoint iteration {} <<<<\n".format(self.FSIIter), self.mpiComm)

            if self.manager.mechanical:
                # --- Solid to fluid mechanical transfer --- #
                self.solidToFluidMechaTransfer()
                self.solidToFluidAdjointTransfer()
                # --- Fluid mesh morphing --- #
                mpiPrint('\nPerforming mesh deformation...\n', self.mpiComm)
                self.meshDefTimer.start()
                self.FluidSolver.meshUpdate(self.step.timeIter)
                self.meshDefTimer.stop()
                self.meshDefTimer.cumul()
            self.FluidSolver.boundaryConditionsUpdate()

            # --- Fluid solver call for FSI subiteration --- #
            mpiPrint('\nLaunching adjoint fluid solver...', self.mpiComm)
            self.fluidSolverTimer.start()
            self.FluidSolver.run(*self.step.timeFrame())
            self.fluidSolverTimer.stop()
            self.fluidSolverTimer.cumul()
            mpiBarrier(self.mpiComm)

            if self.manager.mechanical:
                # --- Fluid to solid mechanical transfer --- #
                mpiPrint('\nProcessing interface fluid adjoint loads...\n', self.mpiComm)
                self.fluidToSolidMechaTransfer()
                self.fluidToSolidAdjointTransfer()
            mpiBarrier(self.mpiComm)

            # --- Solid solver call for FSI subiteration --- #
            mpiPrint('\nLaunching adjoint solid solver...\n', self.mpiComm)
            if self.myid in self.manager.getSolidSolverProcessors():
                self.solidSolverTimer.start()
                self.SolidSolver.run(*self.step.timeFrame())
                self.solidSolverTimer.stop()
                self.solidSolverTimer.cumul()
            self.solidHasRun = True

            if self.manager.mechanical:
                # --- Compute the mechanical residual --- #
                res = self.computeSolidInterfaceAdjointResidual()
                self.errValue = self.criterion.update(res)
                mpiPrint('\nFSI error value : {}\n'.format(self.errValue), self.mpiComm)
            else:
                self.errValue = 0.0
            
            self.errValue_CHT = 0.0

            # --- Monitor the coupling convergence --- #
            self.FSIConv = self.criterion.isVerified(self.errValue, self.errValue_CHT)

            if self.manager.mechanical:
                # --- Relaxe the solid position --- #
                mpiPrint('\nProcessing interface displacements...\n', self.mpiComm)
                self.relaxSolidAdjointLoad()

            # --- Update the solvers for the next BGS steady iteration --- #
            if self.manager.computationType == 'steady':

                if self.myid in self.manager.getSolidSolverProcessors():
                    self.SolidSolver.steadyUpdate()
                self.FluidSolver.steadyUpdate()

            # --- Update the FSI iteration and history --- #
            if self.writeInFSIloop == True:
                self.writeRealTimeData()
            self.FSIIter += 1

        mpiPrint("BGS is Converged",self.mpiComm,titlePrint)

    def fluidToSolidAdjointTransfer(self):

        self.communicationTimer.start()
        self.interfaceInterpolator.getAdjointDisplacementFromFluidSolver()
        self.interfaceInterpolator.interpolateFluidAdjointDisplacementOnSolidMesh()
        self.interfaceInterpolator.setAdjointDisplacementToSolidSolver(self.step.dt)
        self.communicationTimer.stop()
        self.communicationTimer.cumul()

    def solidToFluidAdjointTransfer(self):

        self.communicationTimer.start()
        self.interfaceInterpolator.getAdjointLoadsFromSolidSolver()
        self.interfaceInterpolator.interpolateSolidAdjointLoadsOnFluidMesh()
        self.interfaceInterpolator.setAdjointLoadsToFluidSolver(self.step.dt)
        self.communicationTimer.stop()
        self.communicationTimer.cumul()

    def __unsteadyRun(self):
        RuntimeError("Unsteady adjoint not implemented")

    def computeSolidInterfaceAdjointResidual(self):

        ns = self.interfaceInterpolator.getNs()
        d = self.interfaceInterpolator.getd()

        # --- Get the predicted (computed) solid interface adjoint loads from the solid solver --- #
        predictedAdjointLoad = FlexInterfaceData(ns+d, 3, self.mpiComm)

        if self.myid in self.manager.getSolidInterfaceProcessors():
            localSolidInterfaceAdjointLoad_X, localSolidInterfaceAdjointLoad_Y, localSolidInterfaceAdjointLoad_Z = self.SolidSolver.getNodalAdjointLoads()
            for iVertex in range(self.manager.getNumberOfLocalSolidInterfaceNodes()):
                iGlobalVertex = self.manager.getGlobalIndex('solid', self.myid, iVertex)
                predictedAdjointLoad[iGlobalVertex] = [localSolidInterfaceAdjointLoad_X[iVertex], localSolidInterfaceAdjointLoad_Y[iVertex], localSolidInterfaceAdjointLoad_Z[iVertex]]

        predictedAdjointLoad.assemble()

        # --- Calculate the residual (vector and norm) --- #
        mpiPrint("\nCompute FSI residual based on solid adjoint load displacement.", self.mpiComm)
        self.solidInterfaceResidual.set(predictedAdjointLoad - self.interfaceInterpolator.solidInterfaceAdjointLoads)

        return self.solidInterfaceResidual

    def relaxSolidAdjointLoad(self):

        # --- Set the relaxation parameter --- #
        self.setOmegaMecha()

        # --- Relax the solid interface position --- #
        self.interfaceInterpolator.solidInterfaceAdjointLoads += (self.omegaMecha*self.solidInterfaceResidual)
        