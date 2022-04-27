#! /usr/bin/env python
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
import scipy as sp
import sys
import traceback
import copy
import sys

import ccupydo
from .utilities import *
from .interfaceData import FlexInterfaceData

np.set_printoptions(threshold=sys.maxsize)

# ----------------------------------------------------------------------
#    Algorithm class
# ----------------------------------------------------------------------

class Algorithm(object):
    """
    Des.
    """

    def __init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, deltaT, totTime, timeIterTreshold=-1, dtSave=0, mpiComm=None):
        """
        Des.
        """

        mpiPrint('\n***************************** Initializing FSI algorithm *****************************', mpiComm)

        self.mpiComm = mpiComm
        self.manager = Manager
        self.FluidSolver = FluidSolver
        self.SolidSolver = SolidSolver
        self.interfaceInterpolator = InterfaceInterpolator
        

        self.globalTimer = Timer()
        self.communicationTimer = Timer()
        self.meshDefTimer = Timer()
        self.fluidSolverTimer = Timer()
        self.solidSolverTimer = Timer()
        self.solidRemeshingTimer = Timer()
        self.fluidRemeshingTimer = Timer()

        
        self.deltaT = deltaT
        self.totTime = totTime
        self.dtSave = dtSave
        self.timeIterTreshold = timeIterTreshold

        
        self.time = 0.0
        self.timeIter = 0
        
        if self.mpiComm != None:
            self.myid = self.mpiComm.Get_rank()
            self.mpiSize = self.mpiComm.Get_size()
        else:
            self.myid = 0
            self.mpiSize = 1

        self.solidHasRun = False

    def setFSIInitialConditions(self):
        """
        Des.
        """

        if self.manager.mechanical:
            if self.manager.computationType == 'unsteady':
                if self.myid in self.manager.getSolidSolverProcessors():
                    self.SolidSolver.setInitialDisplacements()
                self.interfaceInterpolator.getDisplacementFromSolidSolver()
                self.interfaceInterpolator.interpolateSolidDisplacementOnFluidMesh()
                self.interfaceInterpolator.setDisplacementToFluidSolver(self.time)
                self.FluidSolver.setInitialMeshDeformation()
            else:
                if self.myid in self.manager.getSolidSolverProcessors():
                    self.SolidSolver.setInitialDisplacements()
                self.interfaceInterpolator.getDisplacementFromSolidSolver()
                self.interfaceInterpolator.interpolateSolidDisplacementOnFluidMesh()
                self.interfaceInterpolator.setDisplacementToFluidSolver(self.time)
                self.FluidSolver.setInitialMeshDeformation()

        if self.manager.thermal:
            if self.interfaceInterpolator.chtTransferMethod == 'hFFB' or self.interfaceInterpolator.chtTransferMethod == 'TFFB':
                self.FluidSolver.setInitialInterfaceHeatFlux()
            elif self.interfaceInterpolator.chtTransferMethod == 'hFTB' or self.interfaceInterpolator.chtTransferMethod == 'FFTB':
                self.FluidSolver.setInitialInterfaceTemperature()
            self.FluidSolver.boundaryConditionsUpdate()

    def fluidToSolidMechaTransfer(self):
        """
        Des.
        """

        self.communicationTimer.start()
        self.interfaceInterpolator.getLoadsFromFluidSolver()
        self.interfaceInterpolator.interpolateFluidLoadsOnSolidMesh()
        self.interfaceInterpolator.setLoadsToSolidSolver(self.time)
        self.communicationTimer.stop()
        self.communicationTimer.cumul()

    def solidToFluidMechaTransfer(self):
        """
        Des.
        """

        self.communicationTimer.start()
        self.interfaceInterpolator.interpolateSolidDisplacementOnFluidMesh()
        self.interfaceInterpolator.setDisplacementToFluidSolver(self.time)
        self.communicationTimer.stop()
        self.communicationTimer.cumul()

    def solidToFluidThermalTransfer(self):
        """
        Des.
        """

        self.communicationTimer.start()
        if self.interfaceInterpolator.chtTransferMethod == 'TFFB' or self.interfaceInterpolator.chtTransferMethod == 'hFFB':
            self.interfaceInterpolator.interpolateSolidHeatFluxOnFluidMesh()
            self.interfaceInterpolator.setHeatFluxToFluidSolver(self.time)
        elif self.interfaceInterpolator.chtTransferMethod == 'FFTB' or self.interfaceInterpolator.chtTransferMethod == 'hFTB':
            self.interfaceInterpolator.interpolateSolidTemperatureOnFluidMesh()
            self.interfaceInterpolator.setTemperatureToFluidSolver(self.time)
        self.communicationTimer.stop()
        self.communicationTimer.cumul()

    def fluidToSolidThermalTransfer(self):
        """
        Des.
        """

        self.communicationTimer.start()
        if self.interfaceInterpolator.chtTransferMethod == 'TFFB':
            self.interfaceInterpolator.getTemperatureFromFluidSolver()
            self.interfaceInterpolator.interpolateFluidTemperatureOnSolidMesh()
            self.interfaceInterpolator.setTemperatureToSolidSolver(self.time)
        elif self.interfaceInterpolator.chtTransferMethod == 'FFTB':
            self.interfaceInterpolator.getHeatFluxFromFluidSolver()
            self.interfaceInterpolator.interpolateFluidHeatFluxOnSolidMesh()
            self.interfaceInterpolator.setHeatFluxToSolidSolver(self.time)
        elif self.interfaceInterpolator.chtTransferMethod == 'hFFB' or self.interfaceInterpolator.chtTransferMethod == 'hFTB':
            self.interfaceInterpolator.getRobinTemperatureFromFluidSolver()
            self.interfaceInterpolator.interpolateFluidRobinTemperatureOnSolidMesh()
            self.interfaceInterpolator.setRobinHeatFluxToSolidSolver(self.time)
        self.communicationTimer.stop()
        self.communicationTimer.cumul()

class AlgorithmExplicit(Algorithm):
    """
    Des.
    """

    def __init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, deltaT, totTime, timeIterTreshold=-1, dtSave=0, mpiComm=None):
        """
        Des.
        """

        Algorithm.__init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, deltaT, totTime, timeIterTreshold, dtSave, mpiComm)

    def run(self):
        """
        Des.
        """
        
        # --- Initialize output manager --- #
        self.iniRealTimeData()

        mpiPrint('\n**********************************', self.mpiComm)
        mpiPrint('*         Begin FSI computation            *', self.mpiComm)
        mpiPrint('**********************************\n', self.mpiComm)

        self.globalTimer.start()

        #If no restart
        mpiPrint('Setting FSI initial conditions...', self.mpiComm)
        self.setFSIInitialConditions()
        mpiPrint('\nFSI initial conditions are set', self.mpiComm)
        
        try:
            if self.manager.computationType == 'unsteady':
                self.__unsteadyRun()
            else:
                raise Exception('Explicit coupling is only valid for unsteady computations!')
        except:
            mpiPrint('\nA DIVINE ERROR OCCURED...EXITING COMPUTATION\n', self.mpiComm)
            traceback.print_exc()
        finally:
            self.globalTimer.stop()
            self.globalTimer.cumul()
            
            mpiBarrier(self.mpiComm)
            
            mpiPrint('\n*************************', self.mpiComm)
            mpiPrint('*    End FSI computation    *', self.mpiComm)
            mpiPrint('*************************\n', self.mpiComm)
            
            self.printExitInfo()
            
            # --- Exit the solid solver --- #
            if self.myid in self.manager.getSolidSolverProcessors():
                self.SolidSolver.exit()

            # --- Exit the fluid solver --- #
            self.FluidSolver.exit()
    
            # --- Exit computation --- #
            mpiBarrier(self.mpiComm)

    def __unsteadyRun(self):
        """
        Des.
        """

        #If no restart
        nbTimeIter = int(self.totTime/self.deltaT-1)
        prevWrite = self.time

        mpiPrint('Begin time integration\n', self.mpiComm)

        # --- External temporal loop --- #
        while self.timeIter <= nbTimeIter:
            
            mpiPrint("\n>>>> Time iteration {} <<<<".format(self.timeIter), self.mpiComm)

            # --- Preprocess the temporal iteration --- #
            self.FluidSolver.preprocessTimeIter(self.timeIter)
            if self.myid in self.manager.getSolidSolverProcessors():
                self.SolidSolver.preprocessTimeIter(self.timeIter)

            # --- Internal FSI loop --- #
            self.fsiCoupling()
            # --- End of FSI loop --- #

            mpiBarrier(self.mpiComm)

            # --- Update the fluid and solid solver for the next time step --- #
            if self.myid in self.manager.getSolidSolverProcessors():
                self.SolidSolver.update()
            self.FluidSolver.update(self.deltaT)

            # --- Write fluid and solid solution, update FSI history  ---#
            if (self.time > 0) and (self.time-prevWrite > self.dtSave):

                prevWrite = self.time
                self.FluidSolver.save(self.timeIter)
                if self.myid in self.manager.getSolidSolverProcessors():
                    self.SolidSolver.save()

            self.writeRealTimeData()

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
            # ---

            self.timeIter += 1
            self.time += self.deltaT
        # --- End of the temporal loop --- #

    def fsiCoupling(self):
        """
        FSI explicit coupling
        """

        if self.timeIter > self.timeIterTreshold:
            mpiPrint('\n*************** Enter explicit FSI coupling ***************', self.mpiComm)

        if self.manager.mechanical:
            # --- Solid to fluid mechanical transfer --- #
            self.interfaceInterpolator.getDisplacementFromSolidSolver()
            self.solidToFluidMechaTransfer()
            # --- Fluid mesh morphing --- #
            mpiPrint('\nPerforming mesh deformation...\n', self.mpiComm)
            self.meshDefTimer.start()
            self.FluidSolver.meshUpdate(self.timeIter)
            self.meshDefTimer.stop()
            self.meshDefTimer.cumul()

        if self.manager.thermal and self.solidHasRun:
            # --- Solid to fluid thermal transfer --- #
            self.solidToFluidThermalTransfer()
        self.FluidSolver.boundaryConditionsUpdate()

        # --- Fluid solver call  --- #
        mpiPrint('\nLaunching fluid solver...', self.mpiComm)
        self.fluidSolverTimer.start()
        self.FluidSolver.run(self.time-self.deltaT, self.time)
        self.fluidSolverTimer.stop()
        self.fluidSolverTimer.cumul()
        mpiBarrier(self.mpiComm)

        if self.timeIter > self.timeIterTreshold:
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
                self.SolidSolver.run(self.time-self.deltaT, self.time)
                self.solidSolverTimer.stop()
                self.solidSolverTimer.cumul()
            self.solidHasRun = True

    def iniRealTimeData(self):
        """
        Des
        """

        if self.myid in self.manager.getSolidSolverProcessors():
            self.SolidSolver.initRealTimeData()

    def writeRealTimeData(self):
        """
        Des
        """

        if self.myid == 0:
            self.FluidSolver.saveRealTimeData(self.time, 0)
            if self.timeIter >= self.timeIterTreshold:
                self.SolidSolver.saveRealTimeData(self.time, 0)

    def printExitInfo(self):
        """
        Des
        """

        mpiPrint('[cpu FSI total]: ' + str(self.globalTimer.cumulTime) + ' s', self.mpiComm)
        mpiPrint('[cpu FSI fluid mesh mapping]: ' + str(self.interfaceInterpolator.mappingTimer.cumulTime) + ' s', self.mpiComm)
        mpiPrint('[cpu FSI fluid mesh deformation]: ' + str(self.meshDefTimer.cumulTime) + ' s', self.mpiComm)
        mpiPrint('[cpu FSI communications]: ' + str(self.communicationTimer.cumulTime) + ' s', self.mpiComm)
        mpiPrint('[cpu FSI fluid solver]: ' + str(self.fluidSolverTimer.cumulTime) + ' s', self.mpiComm)
        mpiPrint('[cpu FSI solid solver]: ' + str(self.solidSolverTimer.cumulTime) + ' s', self.mpiComm)
        mpiPrint('[cpu FSI fluid remeshing]: ' + str(self.fluidRemeshingTimer.cumulTime) + ' s', self.mpiComm)
        mpiPrint('[cpu FSI solid remeshing]: ' + str(self.solidRemeshingTimer.cumulTime) + ' s', self.mpiComm)
        mpiPrint('[Time steps FSI]: ' + str(self.timeIter), self.mpiComm)
        mpiPrint('[Successful Run FSI]: ' + str(self.time >= (self.totTime - 2*self.deltaT)), self.mpiComm) # NB: self.totTime - 2*self.deltaT is the extreme case that can be encountered due to rounding effects!
        mpiPrint('[Mean n. of FSI Iterations]: ' + str(1), self.mpiComm)

        if self.myid == 0 :
            self.FluidSolver.printRealTimeData(self.time, 0)
            self.SolidSolver.printRealTimeData(self.time, 0)

        mpiPrint('RES-FSI-FSIhistory: ' + str(self.timeIter) + '\t' + str(self.time) + '\t' + str(1.0) + '\t' + str(1) + '\n', self.mpiComm)

class AlgorithmBGSStaticRelax(Algorithm):
    """
    Des.
    """

    def __init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, Criterion, nbFSIIterMax, deltaT, totTime, timeIterTreshold=-1, dtSave=0, omegaBoundList=[1.0,1.0], mpiComm=None):
        """
        Des.
        """

        Algorithm.__init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, deltaT, totTime, timeIterTreshold, dtSave, mpiComm)

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
        """
        Des.
        """

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
        """
        Des.
        """

        # --- Initialize interface data --- #
        self.initInterfaceData()
        
        # --- Initialize output manager --- #
        self.iniRealTimeData()

        mpiPrint('\n**********************************', self.mpiComm)
        mpiPrint('*         Begin FSI computation            *', self.mpiComm)
        mpiPrint('**********************************\n', self.mpiComm)

        self.globalTimer.start()

        #If no restart
        mpiPrint('Setting FSI initial conditions...', self.mpiComm)
        self.setFSIInitialConditions()
        mpiPrint('\nFSI initial conditions are set', self.mpiComm)
        
        try:
            if self.manager.computationType == 'unsteady':
                self.__unsteadyRun()
            else:
                self.time = self.totTime
                self.timeIter = 1
                self.deltaT = self.totTime
                self.writeInFSIloop = True
                self.fsiCoupling()
                self.FluidSolver.save(self.timeIter)
                if self.myid in self.manager.getSolidSolverProcessors():
                    self.SolidSolver.save()
                self.totNbOfFSIIt = self.FSIIter
        except:
            mpiPrint('\nA DIVINE ERROR OCCURED...EXITING COMPUTATION\n', self.mpiComm)
            traceback.print_exc()
        finally:
            self.globalTimer.stop()
            self.globalTimer.cumul()
            
            mpiBarrier(self.mpiComm)
            
            mpiPrint('\n*************************', self.mpiComm)
            mpiPrint('*    End FSI computation    *', self.mpiComm)
            mpiPrint('*************************\n', self.mpiComm)
            
            self.printExitInfo()
            
            # --- Exit the solid solver --- #
            if self.myid in self.manager.getSolidSolverProcessors():
                self.SolidSolver.exit()

            # --- Exit the fluid solver --- #
            self.FluidSolver.exit()
    
            # --- Exit computation --- #
            mpiBarrier(self.mpiComm)

    def __unsteadyRun(self):
        """
        Des.
        """

        #If no restart
        nbTimeIter = int(self.totTime/self.deltaT-1)
        prevWrite = self.time

        mpiPrint('Begin time integration\n', self.mpiComm)

        # --- External temporal loop --- #
        while self.timeIter <= nbTimeIter:
            
            mpiPrint("\n>>>> Time iteration {} <<<<".format(self.timeIter), self.mpiComm)

            # --- Preprocess the temporal iteration --- #
            self.FluidSolver.preprocessTimeIter(self.timeIter)
            if self.myid in self.manager.getSolidSolverProcessors():
                self.SolidSolver.preprocessTimeIter(self.timeIter)

            # --- Internal FSI loop --- #
            self.fsiCoupling()
            # --- End of FSI loop --- #

            mpiBarrier(self.mpiComm)
            
            if self.timeIter > 0:
                self.totNbOfFSIIt += self.FSIIter

            # --- Update the fluid and solid solver for the next time step --- #
            if self.myid in self.manager.getSolidSolverProcessors():
                self.SolidSolver.update()
            self.FluidSolver.update(self.deltaT)

            # --- Write fluid and solid solution, update FSI history  ---#

            if (self.time > 0) and (self.time-prevWrite > self.dtSave):

                prevWrite = self.time
                self.FluidSolver.save(self.timeIter)
                if self.myid in self.manager.getSolidSolverProcessors():
                    self.SolidSolver.save()

            self.writeRealTimeData()

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
            # ---

            if self.timeIter >= self.timeIterTreshold and self.predictor:
                # --- Displacement predictor for the next time step and update of the solid solution --- #
                mpiPrint('\nSolid displacement prediction for next time step', self.mpiComm)
                self.solidDisplacementPredictor()

            self.timeIter += 1
            self.time += self.deltaT
        # --- End of the temporal loop --- #

    def iniRealTimeData(self):
        """
        Des
        """

        if self.myid == 0:
            self.FluidSolver.initRealTimeData()
        if self.myid in self.manager.getSolidSolverProcessors():
            self.SolidSolver.initRealTimeData()
        histFile = open('FSIhistory.ascii', "w")
        histFile.write("{0:>12s}   {1:>12s}   {2:>12s}   {3:>12s}   {4:>12s}   {5:>12s}   {6:>12s}\n".format("TimeIter", "Time", "FSIError", "CHTError", "FSINbIter", "omegaMecha", "omegaThermal"))
        histFile.close()

    def writeRealTimeData(self):
        """
        Des
        """

        if self.myid == 0:
            self.FluidSolver.saveRealTimeData(self.time, self.FSIIter)
            if self.timeIter >= self.timeIterTreshold:
                self.SolidSolver.saveRealTimeData(self.time, self.FSIIter)
            histFile = open('FSIhistory.ascii', "a")
            histFile.write('{0:12d}   {1:.6e}   {2:.6e}   {3:.6e}   {4:12d}   {5:.6e}   {6:.6e}\n'.format(self.timeIter, self.time, self.errValue, self.errValue_CHT, self.FSIIter, self.omegaMecha, self.omegaThermal))
            histFile.close()

    def getMeanNbOfFSIIt(self):
        """
        Des
        """
        if self.manager.computationType == 'unsteady':
            if self.timeIter > 1:
                return float(self.totNbOfFSIIt)/(self.timeIter-1)
            else:
                return 0.0
        else:
            return self.FSIIter

    def printExitInfo(self):
        """
        Des
        """

        mpiPrint('[cpu FSI total]: ' + str(self.globalTimer.cumulTime) + ' s', self.mpiComm)
        mpiPrint('[cpu FSI fluid mesh mapping]: ' + str(self.interfaceInterpolator.mappingTimer.cumulTime) + ' s', self.mpiComm)
        mpiPrint('[cpu FSI fluid mesh deformation]: ' + str(self.meshDefTimer.cumulTime) + ' s', self.mpiComm)
        mpiPrint('[cpu FSI communications]: ' + str(self.communicationTimer.cumulTime) + ' s', self.mpiComm)
        mpiPrint('[cpu FSI fluid solver]: ' + str(self.fluidSolverTimer.cumulTime) + ' s', self.mpiComm)
        mpiPrint('[cpu FSI solid solver]: ' + str(self.solidSolverTimer.cumulTime) + ' s', self.mpiComm)
        mpiPrint('[cpu FSI fluid remeshing]: ' + str(self.fluidRemeshingTimer.cumulTime) + ' s', self.mpiComm)
        mpiPrint('[cpu FSI solid remeshing]: ' + str(self.solidRemeshingTimer.cumulTime) + ' s', self.mpiComm)
        mpiPrint('[Time steps FSI]: ' + str(self.timeIter), self.mpiComm)
        mpiPrint('[Successful Run FSI]: ' + str(self.time >= (self.totTime - 2*self.deltaT)), self.mpiComm) # NB: self.totTime - 2*self.deltaT is the extreme case that can be encountered due to rounding effects!
        mpiPrint('[Mean n. of FSI Iterations]: ' + str(self.getMeanNbOfFSIIt()), self.mpiComm)

        if self.myid == 0 :
            self.FluidSolver.printRealTimeData(self.time, self.FSIIter)
            self.SolidSolver.printRealTimeData(self.time, self.FSIIter)

        mpiPrint('RES-FSI-FSIhistory: ' + str(self.timeIter) + '\t' + str(self.time) + '\t' + str(self.errValue) + '\t' + str(self.FSIIter) + '\n', self.mpiComm)

    def fsiCoupling(self):
        """
        Block Gauss Seidel (BGS) method for strong coupling FSI
        """

        if self.timeIter > self.timeIterTreshold:
            nbFSIIter = self.nbFSIIterMax
            mpiPrint('\n*************** Enter Block Gauss Seidel (BGS) method for strong coupling FSI ***************', self.mpiComm)
        else:
             nbFSIIter = 1

        self.FSIIter = 0
        self.FSIConv = False
        self.errValue = 1e12
        self.errValue_CHT = 1e6

        while ((self.FSIIter < nbFSIIter) and (not self.criterion.isVerified(self.errValue, self.errValue_CHT))):
            mpiPrint("\n>>>> FSI iteration {} <<<<\n".format(self.FSIIter), self.mpiComm)

            if self.manager.mechanical:
                # --- Solid to fluid mechanical transfer --- #
                self.solidToFluidMechaTransfer()
                # --- Fluid mesh morphing --- #
                mpiPrint('\nPerforming mesh deformation...\n', self.mpiComm)
                self.meshDefTimer.start()
                self.FluidSolver.meshUpdate(self.timeIter)
                self.meshDefTimer.stop()
                self.meshDefTimer.cumul()
            if self.manager.thermal and self.solidHasRun:
                # --- Solid to fluid thermal transfer --- #
                self.solidToFluidThermalTransfer()
            self.FluidSolver.boundaryConditionsUpdate()

            # --- Fluid solver call for FSI subiteration --- #
            mpiPrint('\nLaunching fluid solver...', self.mpiComm)
            self.fluidSolverTimer.start()
            self.FluidSolver.run(self.time-self.deltaT, self.time)
            self.fluidSolverTimer.stop()
            self.fluidSolverTimer.cumul()
            mpiBarrier(self.mpiComm)

            if self.timeIter > self.timeIterTreshold:
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
                    self.SolidSolver.run(self.time-self.deltaT, self.time)
                    self.solidSolverTimer.stop()
                    self.solidSolverTimer.cumul()
                self.solidHasRun = True

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

                # --- Monitor the coupling convergence --- #
                self.FSIConv = self.criterion.isVerified(self.errValue, self.errValue_CHT)

                if self.manager.mechanical:
                    # --- Relaxe the solid position --- #
                    mpiPrint('\nProcessing interface displacements...\n', self.mpiComm)
                    self.relaxSolidPosition()

                if self.manager.thermal:
                    # --- Relaxe thermal data --- #
                    self.relaxCHT()

            if self.writeInFSIloop == True:
                self.writeRealTimeData()

            self.FSIIter += 1
            if self.manager.computationType != 'unsteady':
                self.time += self.deltaT

            # --- Update the solvers for the next BGS iteration --- #
            if self.myid in self.manager.getSolidSolverProcessors():
                self.SolidSolver.bgsUpdate()
            self.FluidSolver.bgsUpdate()

        if self.timeIter > self.timeIterTreshold:
            mpiPrint('\n*************** BGS is converged ***************', self.mpiComm)

    def computeSolidInterfaceResidual(self):
        """
        Des.
        """

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
        """
        Des.
        """

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
        """
        Des
        """

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

        # --- Predict the solid position for the next time step --- #
        if self.predictorOrder == 1:
            mpiPrint("First order predictor.", self.mpiComm)
            self.interfaceInterpolator.solidInterfaceDisplacement += (self.alpha_0*self.deltaT*self.solidInterfaceVelocity)
        else:
            mpiPrint("Second order predictor.", self.mpiComm)
            self.interfaceInterpolator.solidInterfaceDisplacement += (self.alpha_0*self.deltaT*self.solidInterfaceVelocity + self.alpha_1*self.deltaT*(self.solidInterfaceVelocity-self.solidInterfaceVelocitynM1))

    def setOmegaMecha(self):
        """
        Des.
        """

        self.omegaMecha = self.omegaBoundMecha

        mpiPrint('Static under-relaxation summary, mechanical : {}'.format(self.omegaMecha), self.mpiComm)


    def setOmegaThermal(self):
        """
        Des.
        """

        self.omegaThermal = self.omegaBoundThermal

        mpiPrint('Static under-relaxation summary, thermal : {}'.format(self.omegaThermal), self.mpiComm)

    def relaxSolidPosition(self):
        """
        Des.
        """

        # --- Set the relaxation parameter --- #
        self.setOmegaMecha()

        # --- Relax the solid interface position --- #
        self.interfaceInterpolator.solidInterfaceDisplacement += (self.omegaMecha*self.solidInterfaceResidual)

    def relaxCHT(self):
        """
        Des.
        """

        self.setOmegaThermal()

        if self.interfaceInterpolator.chtTransferMethod == 'hFFB' or self.interfaceInterpolator.chtTransferMethod == 'TFFB':
            self.interfaceInterpolator.solidInterfaceHeatFlux += (self.omegaThermal*self.solidHeatFluxResidual)
        elif self.interfaceInterpolator.chtTransferMethod == 'hFTB' or self.interfaceInterpolator.chtTransferMethod == 'FFTB':
            self.interfaceInterpolator.solidInterfaceTemperature += (self.omegaThermal*self.solidTemperatureResidual)

class AlgorithmBGSAitkenRelax(AlgorithmBGSStaticRelax):

    def __init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, Criterion, nbFSIIterMax, deltaT, totTime, timeIterTreshold=-1, dtSave=0, omegaBoundList=[1.0, 1.0], mpiComm=None):
        """
        Des.
        """

        AlgorithmBGSStaticRelax.__init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, Criterion, nbFSIIterMax, deltaT, totTime, timeIterTreshold, dtSave, omegaBoundList, mpiComm)


        self.solidInterfaceResidualkM1 = None
        self.solidHeatFluxResidualkM1 = None
        self.solidTemperatureResidualkM1 = None
        self.aitkenCritMecha = 'max'
        self.aitkenCritThermal = 'max'

    def initInterfaceData(self):
        """
        Des.
        """

        AlgorithmBGSStaticRelax.initInterfaceData(self)
        ns = self.interfaceInterpolator.getNs()
        d = self.interfaceInterpolator.getd()

        self.solidInterfaceResidualkM1 = FlexInterfaceData(ns+d, 3, self.mpiComm)
        self.solidHeatFluxResidualkM1 = FlexInterfaceData(ns+d, 3, self.mpiComm)
        self.solidTemperatureResidualkM1 = FlexInterfaceData(ns+d, 1, self.mpiComm)



    def setOmegaMecha(self):
        """
        Des.
        """

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
        """
        Des.
        """

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

class AlgorithmIQN_ILS(AlgorithmBGSAitkenRelax):
    """
    Des.
    """

    def __init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, Criterion, nbFSIIterMax, deltaT, totTime, timeIterTreshold=-1, dtSave=0, omegaBoundList= [1.0, 1.0], nbTimeToKeep=0, computeTangentMatrixBasedOnFirstIt = False, mpiComm=None):
        """
        Des.
        """

        AlgorithmBGSAitkenRelax.__init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, Criterion, nbFSIIterMax, deltaT, totTime, timeIterTreshold, dtSave, omegaBoundList, mpiComm)

        # --- Number of previous time steps used in the approximation of the tangent matrix --- #
        self.nbTimeToKeep = nbTimeToKeep

        # --- Option which allows to build the tangent matrix of a given time step using differences with respect to the first FSI iteration (delta_r_k = r_k+1 - r_0) instead of the previous iteration (delta_r_k = r_k+1 - r_k) --- #
        self.computeTangentMatrixBasedOnFirstIt = computeTangentMatrixBasedOnFirstIt

        # --- Option which determines the way the c coefficients are computed either using Degroote's QR decompoistion or simply using np.linalg.lstsq
        self.useQR = True
        self.tollQR = 1.0e-1 # Tolerance employed for the QR decomposition
        self.qrFilter = 'Haelterman' # Type of QR filtering employed. Possible choices are 'Degroote1', 'Degroote2', and 'Haelterman' (see 'qrSolve()' function below)
        
        self.maxNbOfItReached = False
        self.convergenceReachedInOneIt = False
        
        # --- Global V and W matrices for IQN-ILS algorithm, including information from previous time steps --- #
        self.V = []
        self.W = []
    
    def qrSolve(self, V, W, res):
        
        if self.qrFilter == 'Degroote1': # QR filtering as described by J. Degroote et al. Computers and Structures, 87, 793-801 (2009).
            Q, R = sp.linalg.qr(V, mode='economic')
            s = np.dot(np.transpose(Q), -res)
            toll = self.tollQR*sp.linalg.norm(R, 2)
            c = solve_upper_triangular_mod(R, s, toll)
        
        elif self.qrFilter == 'Degroote2': # QR filtering as described by J. Degroote et al. CMAME, 199, 2085-2098 (2010).
            Q, R, V, W = QRfiltering(V, W, self.tollQR)
            s = np.dot(np.transpose(Q), -res)
            c = np.linalg.solve(R, s)
        
        elif self.qrFilter == 'Haelterman': # 'Modified' QR filtering as described by R. Haelterman et al. Computers and Structures, 171, 9-17 (2016).
            Q, R, V, W = QRfiltering_mod(V, W, self.tollQR)
            s = np.dot(np.transpose(Q), -res)
            c = np.linalg.solve(R, s)
        
        else:
            raise NameError('IQN-ILS Algorithm: the QR filtering technique is unknown!')
        
        return c, W
    
    def fsiCoupling(self):
        """
        Interface Quasi Newton - Inverse Least Square (IQN-ILS) method for strong coupling FSI
        """

        if self.timeIter > self.timeIterTreshold:
            nbFSIIter = self.nbFSIIterMax
            mpiPrint('\n*************** Enter Interface Quasi Newton - Inverse Least Square (IQN-ILS) method for strong coupling FSI ***************', self.mpiComm)
        else:
             nbFSIIter = 1

        self.FSIIter = 0
        self.FSIConv = False
        self.errValue = 1.0
        self.errValue_CHT = 1e6 # Just for compatibility. CHT not yet implemented for the IQN-ILS algorithm.
        
        ns = self.interfaceInterpolator.getNs()
        d = self.interfaceInterpolator.getd()

        # --- Initialize all the quantities used in the IQN-ILS method --- #
        res = FlexInterfaceData(ns, 3, self.mpiComm)
        solidInterfaceResidual0 = FlexInterfaceData(ns+d, 3, self.mpiComm)

        solidInterfaceDisplacement_tilde = FlexInterfaceData(ns, 3, self.mpiComm)
        solidInterfaceDisplacement_tilde1 = FlexInterfaceData(ns, 3, self.mpiComm)

        delta_ds = FlexInterfaceData(ns+d, 3, self.mpiComm)

        Vk_mat = np.zeros((self.manager.nDim*ns,1))
        Wk_mat = np.zeros((self.manager.nDim*ns,1))

        delta_ds_loc_X = np.zeros(0)
        delta_ds_loc_Y = np.zeros(0)
        delta_ds_loc_Z = np.zeros(0)

        if (self.nbTimeToKeep!=0 and self.timeIter > 1): # If information from previous time steps is re-used then Vk = V, Wk = W
            Vk = copy.deepcopy(self.V)
            Wk = copy.deepcopy(self.W)
        else: # If information from previous time steps is not re-used then Vk and Wk are empty lists of np.array()
            Vk = []
            Wk = []
        
        nIt = 0

        while ((self.FSIIter < nbFSIIter) and (not self.criterion.isVerified(self.errValue,self.errValue_CHT))):
            mpiPrint("\n>>>> FSI iteration {} <<<<\n".format(self.FSIIter), self.mpiComm)

            # --- Solid to fluid mechanical transfer --- #
            self.solidToFluidMechaTransfer()
            # --- Fluid mesh morphing --- #
            mpiPrint('\nPerforming mesh deformation...\n', self.mpiComm)
            self.meshDefTimer.start()
            self.FluidSolver.meshUpdate(self.timeIter)
            self.meshDefTimer.stop()
            self.meshDefTimer.cumul()

            # --- Fluid solver call for FSI subiteration --- #
            mpiPrint('\nLaunching fluid solver...', self.mpiComm)
            self.fluidSolverTimer.start()
            self.FluidSolver.run(self.time-self.deltaT, self.time)
            self.fluidSolverTimer.stop()
            self.fluidSolverTimer.cumul()
            mpiBarrier(self.mpiComm)

            if self.timeIter > self.timeIterTreshold:
                # --- Fluid to solid mechanical transfer --- #
                mpiPrint('\nProcessing interface fluid loads...\n', self.mpiComm)
                self.fluidToSolidMechaTransfer()
                mpiBarrier(self.mpiComm)

                # --- Solid solver call for FSI subiteration --- #
                mpiPrint('\nLaunching solid solver...\n', self.mpiComm)
                if self.myid in self.manager.getSolidSolverProcessors():
                    self.solidSolverTimer.start()
                    self.SolidSolver.run(self.time-self.deltaT, self.time)
                    self.solidSolverTimer.stop()
                    self.solidSolverTimer.cumul()

                # --- Compute and monitor the FSI residual --- #
                res = self.computeSolidInterfaceResidual()
                self.errValue = self.criterion.update(res)
                mpiPrint('\nFSI error value : {}\n'.format(self.errValue), self.mpiComm)
                self.FSIConv = self.criterion.isVerified(self.errValue)

                # --- Initialize d_tilde for the construction of the Wk matrix -- #
                if self.myid in self.manager.getSolidInterfaceProcessors():
                    localSolidInterfaceDisp_X, localSolidInterfaceDisp_Y, localSolidInterfaceDisp_Z = self.SolidSolver.getNodalDisplacements()
                    for iVertex in range(self.manager.getNumberOfLocalSolidInterfaceNodes()):
                        iGlobalVertex = self.manager.getGlobalIndex('solid', self.myid, iVertex)
                        solidInterfaceDisplacement_tilde[iGlobalVertex] = [localSolidInterfaceDisp_X[iVertex], localSolidInterfaceDisp_Y[iVertex], localSolidInterfaceDisp_Z[iVertex]]

                solidInterfaceDisplacement_tilde.assemble()
                
                if ((self.FSIIter == 0 and (self.nbTimeToKeep == 0 or (self.nbTimeToKeep != 0 and (self.maxNbOfItReached or self.convergenceReachedInOneIt or self.timeIter == 1)))) or self.timeIter < 1): # If information from previous time steps is re-used then this step is only performed at the first iteration of the first time step, otherwise it is performed at the first iteration of every time step
                    # --- Relax the solid position --- #
                    mpiPrint('\nProcessing interface displacements...\n', self.mpiComm)
                    self.relaxSolidPosition()
                else:
                    # --- Construct Vk and Wk matrices for the computation of the approximated tangent matrix --- #
                    mpiPrint('\nCorrect solid interface displacements using IQN-ILS method...\n', self.mpiComm)
                    
                    # --- Start gathering on root process --- #
                    res_X_Gat, res_Y_Gat, res_Z_Gat = mpiGatherInterfaceData(res, ns+d, self.mpiComm, 0)
                    solidInterfaceResidual0_X_Gat, solidInterfaceResidual0_Y_Gat, solidInterfaceResidual0_Z_Gat = mpiGatherInterfaceData(solidInterfaceResidual0, ns+d, self.mpiComm, 0)
                    solidInterfaceDisplacement_tilde_X_Gat, solidInterfaceDisplacement_tilde_Y_Gat, solidInterfaceDisplacement_tilde_Z_Gat = mpiGatherInterfaceData(solidInterfaceDisplacement_tilde, ns, self.mpiComm, 0)
                    solidInterfaceDisplacement_tilde1_X_Gat, solidInterfaceDisplacement_tilde1_Y_Gat, solidInterfaceDisplacement_tilde1_Z_Gat = mpiGatherInterfaceData(solidInterfaceDisplacement_tilde1, ns, self.mpiComm, 0)
                    if self.myid == 0:
                        res_X_Gat_C = res_X_Gat[:ns] # Copies for operating on, length=ns, not ns+d
                        res_Y_Gat_C = res_Y_Gat[:ns]
                        res_Z_Gat_C = res_Z_Gat[:ns]
                        solidInterfaceResidual0_X_Gat_C = solidInterfaceResidual0_X_Gat[:ns] # Copies for operating on, length=ns, not ns+d
                        solidInterfaceResidual0_Y_Gat_C = solidInterfaceResidual0_Y_Gat[:ns]
                        solidInterfaceResidual0_Z_Gat_C = solidInterfaceResidual0_Z_Gat[:ns]
                        if self.FSIIter > 0: # Either information from previous time steps is re-used or not, Vk and Wk matrices are enriched only starting from the second iteration of every FSI loop
                            if self.manager.nDim == 3:
                                delta_res = np.concatenate([res_X_Gat_C - solidInterfaceResidual0_X_Gat_C, res_Y_Gat_C - solidInterfaceResidual0_Y_Gat_C, res_Z_Gat_C - solidInterfaceResidual0_Z_Gat_C], axis=0)
                                delta_d = np.concatenate([solidInterfaceDisplacement_tilde_X_Gat - solidInterfaceDisplacement_tilde1_X_Gat, solidInterfaceDisplacement_tilde_Y_Gat - solidInterfaceDisplacement_tilde1_Y_Gat, solidInterfaceDisplacement_tilde_Z_Gat - solidInterfaceDisplacement_tilde1_Z_Gat], axis = 0)
                            else:
                                delta_res = np.concatenate([res_X_Gat_C - solidInterfaceResidual0_X_Gat_C, res_Y_Gat_C - solidInterfaceResidual0_Y_Gat_C], axis=0)
                                delta_d = np.concatenate([solidInterfaceDisplacement_tilde_X_Gat - solidInterfaceDisplacement_tilde1_X_Gat, solidInterfaceDisplacement_tilde_Y_Gat - solidInterfaceDisplacement_tilde1_Y_Gat], axis = 0)
                            
                            Vk.insert(0, delta_res)
                            Wk.insert(0, delta_d)
                            
                            nIt+=1
                        
                        Vk_mat = np.vstack(Vk).T
                        Wk_mat = np.vstack(Wk).T
                        
                        if (Vk_mat.shape[1] > self.manager.nDim*ns and self.qrFilter == 'Degroote1'): # Remove extra columns if number of iterations (i.e. columns of Vk and Wk) is larger than number of interface degrees of freedom 
                            mpiPrint('WARNING: IQN-ILS Algorithm using \'Degroote1\' QR filter. Approximated stiffness matrix number of columns exceeds the number of degrees of freedom at FSI interface. Extra columns (the oldest ones!) are deleted for next iterations to avoid overdetermined problem!', self.mpiComm)
                            Vk_mat = np.delete(Vk_mat, np.s_[(self.manager.nDim*ns-Vk_mat.shape[1]):], 1)
                            Wk_mat = np.delete(Wk_mat, np.s_[(self.manager.nDim*ns-Wk_mat.shape[1]):], 1)
                        
                        dummy_V = Vk_mat.copy()
                        dummy_W = Wk_mat.copy()
                        
                        if self.manager.nDim == 3:
                            dummy_Res = np.concatenate([res_X_Gat_C, res_Y_Gat_C, res_Z_Gat_C], axis=0)
                        else:
                            dummy_Res = np.concatenate([res_X_Gat_C, res_Y_Gat_C], axis=0)
                        
                        if self.useQR: # Technique described by Degroote et al.
                            c, dummy_W = self.qrSolve(dummy_V, dummy_W, dummy_Res)
                        else:
                            c = np.linalg.lstsq(dummy_V, -dummy_Res)[0] # Classical QR decomposition: NOT RECOMMENDED!
                        
                        if self.manager.nDim == 3:
                            delta_ds_loc = np.split((np.dot(dummy_W,c).T + np.concatenate([res_X_Gat_C, res_Y_Gat_C, res_Z_Gat_C], axis=0)),3,axis=0)
                            
                            delta_ds_loc_X = delta_ds_loc[0]
                            delta_ds_loc_Y = delta_ds_loc[1]
                            delta_ds_loc_Z = delta_ds_loc[2]
                        else:
                            delta_ds_loc = np.split((np.dot(dummy_W,c).T + np.concatenate([res_X_Gat_C, res_Y_Gat_C], axis=0)),2,axis=0)
                            
                            delta_ds_loc_X = delta_ds_loc[0]
                            delta_ds_loc_Y = delta_ds_loc[1]
                            delta_ds_loc_Z = np.zeros(ns)
                        
                        for iVertex in range(delta_ds_loc_X.shape[0]):
                            iGlobalVertex = self.manager.getGlobalIndex('solid', self.myid, iVertex)
                            delta_ds[iGlobalVertex] = [delta_ds_loc_X[iVertex], delta_ds_loc_Y[iVertex], delta_ds_loc_Z[iVertex]]
                    
                    # --- Go back to parallel run --- #
                    mpiBarrier(self.mpiComm)
                    delta_ds.assemble()
                    self.interfaceInterpolator.solidInterfaceDisplacement += delta_ds
                
                if self.computeTangentMatrixBasedOnFirstIt:
                    if self.FSIIter == 0:
                        res.copy(solidInterfaceResidual0)
                        solidInterfaceDisplacement_tilde.copy(solidInterfaceDisplacement_tilde1)
                else:
                    res.copy(solidInterfaceResidual0)
                    solidInterfaceDisplacement_tilde.copy(solidInterfaceDisplacement_tilde1)
            
            if self.writeInFSIloop == True:
                self.writeRealTimeData()
            
            self.FSIIter += 1
        
        # if comm.myself == rootProcess
        
        # update of the matrices V and W at the end of the while
        if self.nbTimeToKeep != 0 and self.timeIter >= 1:
            
            # --- Trick to avoid breaking down of the simulation in the rare cases when, in the initial time steps, FSI convergence is reached without iterating (e.g. starting from a steady condition and using very small time steps), leading to empty V and W matrices ---
            if not (self.FSIIter == 1 and self.FSIConv and len(self.V)==0):
                
                self.convergenceReachedInOneIt = False
                
                # --- Managing situations where FSI convergence is not reached ---
                if (self.FSIIter >= nbFSIIter and not self.FSIConv):
                    mpiPrint('WARNING: IQN-ILS using information from {} previous time steps reached max number of iterations. Next time step is run without using any information from previous time steps!'.format(self.nbTimeToKeep), self.mpiComm)
                    
                    self.maxNbOfItReached = True
                    self.V = []
                    self.W = []
                else:
                    self.maxNbOfItReached = False
                    
                    mpiPrint('\nUpdating V and W matrices...\n', self.mpiComm)
                    
                    self.V.insert(0, Vk_mat[:,0:nIt].T)
                    self.W.insert(0, Wk_mat[:,0:nIt].T)
                    
                    if (self.timeIter > self.nbTimeToKeep and len(self.V) > self.nbTimeToKeep):
                        del self.V[-1]
                        del self.W[-1]
                # --- 
            else:
                mpiPrint('\nWARNING: IQN-ILS algorithm convergence reached in one iteration at the beginning of the simulation. V and W matrices cannot be built. BGS will be employed for the next time step!\n', self.mpiComm)
                self.convergenceReachedInOneIt = True
            # ---
            
        # --- Update the FSI history file --- #
        if self.timeIter > self.timeIterTreshold:
            mpiPrint('\n*************** IQN-ILS is converged ***************', self.mpiComm)

# --- Solid test algorithm ---
class FsiSolidTestAlgorithm(object):
    def __init__(self, _solid):
        self.solid = _solid

    def run(self):
        # --------------------------
        # fake FSI solver
        # --------------------------

        t1 = 0.0  # initial time
        dt = 0.5  # time step size
        nt = 10

        # we want nt time steps
        for j in range(nt):

            # each time step is arbitrarily calculated twice (for testing purpose)
            for i in range(2):

                t2=t1+dt  # time to be calculated

                self.solid.fakeFluidSolver(t2)  # creates some dummy loads for time=t2

                # run solid solver
                print('='*80)
                print("running from %f to %f: try #%d" % (t1,t2,i+1))
                print('='*80)
                self.solid.run(t1,t2)

                # gets the deformed interface
                dx, dy, dz = self.solid.getNodalDisplacements()
                print(dx)
                print(dy)
                print(dz)

            self.solid.update()
            self.solid.save()

            t1=t2 # fsi loop has converged - time t2 is accepted

        # end.
class AlgorithmBGSStaticRelaxAdjoint(AlgorithmBGSStaticRelax):
    def run(self):

        # --- Initialize interface data --- #
        self.initInterfaceData()
        
        # --- Initialize output manager --- #
        self.iniRealTimeData()

        mpiPrint('\n*****************************************', self.mpiComm)
        mpiPrint('*     Begin FSI Adjoint computation     *', self.mpiComm)
        mpiPrint('*****************************************\n', self.mpiComm)

        self.globalTimer.start()

        #If no restart
        mpiPrint('Setting FSI initial conditions...', self.mpiComm)
        self.setFSIInitialConditions()
        mpiPrint('\nFSI initial conditions are set', self.mpiComm)
        
        try:
            self.time = self.totTime
            self.timeIter = 1
            self.deltaT = self.totTime
            self.writeInFSIloop = True
            self.fsiCoupling()
            self.FluidSolver.save(self.timeIter)
            if self.myid in self.manager.getSolidSolverProcessors():
                self.SolidSolver.save()
            self.totNbOfFSIIt = self.FSIIter
        except:
            mpiPrint('\nA DIVINE ERROR OCCURED...EXITING COMPUTATION\n', self.mpiComm)
            traceback.print_exc()
        finally:
            self.globalTimer.stop()
            self.globalTimer.cumul()
            
            mpiBarrier(self.mpiComm)
            
            mpiPrint('\n*****************************************', self.mpiComm)
            mpiPrint('*      End FSI Adjoint computation      *', self.mpiComm)
            mpiPrint('*****************************************\n', self.mpiComm)
            
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

        if self.timeIter > self.timeIterTreshold:
            nbFSIIter = self.nbFSIIterMax
            mpiPrint('\n*************** Enter Block Gauss Seidel (BGS) method for strong coupling FSI ***************', self.mpiComm)
        else:
            nbFSIIter = 1

        self.FSIIter = 0
        self.FSIConv = False
        self.errValue = 1e12
        self.errValue_CHT = 1e6

        while ((self.FSIIter < nbFSIIter) and (not self.criterion.isVerified(self.errValue, self.errValue_CHT))):
            mpiPrint("\n>>>> FSI Adjoint iteration {} <<<<\n".format(self.FSIIter), self.mpiComm)

            if self.manager.mechanical:
                # --- Solid to fluid mechanical transfer --- #
                self.solidToFluidMechaTransfer()
                self.solidToFluidAdjointTransfer()
                # --- Fluid mesh morphing --- #
                mpiPrint('\nPerforming mesh deformation...\n', self.mpiComm)
                self.meshDefTimer.start()
                self.FluidSolver.meshUpdate(self.timeIter)
                self.meshDefTimer.stop()
                self.meshDefTimer.cumul()
            self.FluidSolver.boundaryConditionsUpdate()

            # --- Fluid solver call for FSI subiteration --- #
            mpiPrint('\nLaunching adjoint fluid solver...', self.mpiComm)
            self.fluidSolverTimer.start()
            self.FluidSolver.run(self.time-self.deltaT, self.time)
            self.fluidSolverTimer.stop()
            self.fluidSolverTimer.cumul()
            mpiBarrier(self.mpiComm)

            if self.timeIter > self.timeIterTreshold:
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
                    self.SolidSolver.run(self.time-self.deltaT, self.time)
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

            if self.writeInFSIloop == True:
                self.writeRealTimeData()

            self.FSIIter += 1
            if self.manager.computationType != 'unsteady':
                self.time += self.deltaT

            # --- Update the solvers for the next BGS iteration --- #
            if self.myid in self.manager.getSolidSolverProcessors():
                self.SolidSolver.bgsUpdate()
            self.FluidSolver.bgsUpdate()

        if self.timeIter > self.timeIterTreshold:
            mpiPrint('\n*************** BGS is converged ***************', self.mpiComm)

    def fluidToSolidAdjointTransfer(self):
        """
        Des.
        """

        self.communicationTimer.start()
        self.interfaceInterpolator.getAdjointDisplacementFromFluidSolver()
        self.interfaceInterpolator.interpolateFluidAdjointDisplacementOnSolidMesh()
        self.interfaceInterpolator.setAdjointDisplacementToSolidSolver(self.time)
        self.communicationTimer.stop()
        self.communicationTimer.cumul()

    def solidToFluidAdjointTransfer(self):
        """
        Des.
        """

        self.communicationTimer.start()
        self.interfaceInterpolator.getAdjointLoadsFromSolidSolver()
        self.interfaceInterpolator.interpolateSolidAdjointLoadsOnFluidMesh()
        self.interfaceInterpolator.setAdjointLoadsToFluidSolver(self.time)
        self.communicationTimer.stop()
        self.communicationTimer.cumul()

    def __unsteadyRun(self):
        RuntimeError("Unsteady adjoint not implemented")

    def computeSolidInterfaceAdjointResidual(self):
        """
        Des.
        """

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
        """
        Des.
        """

        # --- Set the relaxation parameter --- #
        self.setOmegaMecha()

        # --- Relax the solid interface position --- #
        self.interfaceInterpolator.solidInterfaceAdjointLoads += (self.omegaMecha*self.solidInterfaceResidual)

class AlgorithmIQN_MVJ(AlgorithmBGSAitkenRelax):
    """
    Des.
    """

    def __init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, Criterion, nbFSIIterMax, deltaT, totTime, timeIterTreshold=-1, dtSave=0, omegaBoundList= [1.0, 1.0], computeTangentMatrixBasedOnFirstIt = False, mpiComm=None):
        """
        Des.
        """

        AlgorithmBGSAitkenRelax.__init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, Criterion, nbFSIIterMax, deltaT, totTime, timeIterTreshold, dtSave, omegaBoundList, mpiComm)

        # --- Indicate if a BGS iteration must be performed --- #
        self.makeBGS = True

        # --- Internal variables for convcergence check --- #
        self.maxNbOfItReached = False
        self.convergenceReachedInOneIt = False

        # --- Option which allows to build the tangent matrix of a given time step using differences with respect to the first FSI iteration (delta_r_k = r_k+1 - r_0) instead of the previous iteration (delta_r_k = r_k+1 - r_k) --- #
        self.computeTangentMatrixBasedOnFirstIt = computeTangentMatrixBasedOnFirstIt
        
        # --- Global V and W matrices for IQN-MVJ algorithm --- #
        self.V = []
        self.W = []

        # --- Curent inverse approximate Jacobian (J^-1 - I) --- #
        ns = self.interfaceInterpolator.getNs()
        self.invJ = np.zeros((self.manager.nDim*ns,self.manager.nDim*ns))
    
    def fsiCoupling(self):
        """
        Interface Quasi Newton - Multi Vector Jacobian (IQN-MVJ) method for strong coupling FSI
        """

        if self.timeIter > self.timeIterTreshold:
            nbFSIIter = self.nbFSIIterMax
            mpiPrint('\n*************** Enter Interface Quasi Newton - Multi Vector Jacobian (IQN-MVJ) method for strong coupling FSI ***************', self.mpiComm)
        else:
             nbFSIIter = 1

        self.FSIIter = 0
        self.FSIConv = False
        self.errValue = 1.0
        self.errValue_CHT = 1e6 # Just for compatibility
        
        ns = self.interfaceInterpolator.getNs()
        d = self.interfaceInterpolator.getd()

        # --- Initialize all the quantities used in the IQN-ILS method --- #
        res = FlexInterfaceData(ns, 3, self.mpiComm)
        solidInterfaceResidual0 = FlexInterfaceData(ns+d, 3, self.mpiComm)

        solidInterfaceDisplacement_tilde = FlexInterfaceData(ns, 3, self.mpiComm)
        solidInterfaceDisplacement_tilde1 = FlexInterfaceData(ns, 3, self.mpiComm)

        delta_ds = FlexInterfaceData(ns+d, 3, self.mpiComm)

        delta_ds_loc_X = np.zeros(0)
        delta_ds_loc_Y = np.zeros(0)
        delta_ds_loc_Z = np.zeros(0)

        # --- Initialises the previous approximate inverse Jacobian and V, W --- #
        Vk = []
        Wk = []
        self.invJprev = self.invJ.copy()
        nIt = 0

        while ((self.FSIIter < nbFSIIter) and (not self.criterion.isVerified(self.errValue,self.errValue_CHT))):
            mpiPrint("\n>>>> FSI iteration {} <<<<\n".format(self.FSIIter), self.mpiComm)

            # --- Solid to fluid mechanical transfer --- #
            self.solidToFluidMechaTransfer()
            # --- Fluid mesh morphing --- #
            mpiPrint('\nPerforming mesh deformation...\n', self.mpiComm)
            self.meshDefTimer.start()
            self.FluidSolver.meshUpdate(self.timeIter)
            self.meshDefTimer.stop()
            self.meshDefTimer.cumul()

            # --- Fluid solver call for FSI subiteration --- #
            mpiPrint('\nLaunching fluid solver...', self.mpiComm)
            self.fluidSolverTimer.start()
            self.FluidSolver.run(self.time-self.deltaT, self.time)
            self.fluidSolverTimer.stop()
            self.fluidSolverTimer.cumul()
            mpiBarrier(self.mpiComm)

            if self.timeIter > self.timeIterTreshold:
                # --- Fluid to solid mechanical transfer --- #
                mpiPrint('\nProcessing interface fluid loads...\n', self.mpiComm)
                self.fluidToSolidMechaTransfer()
                mpiBarrier(self.mpiComm)

                # --- Solid solver call for FSI subiteration --- #
                mpiPrint('\nLaunching solid solver...\n', self.mpiComm)
                if self.myid in self.manager.getSolidSolverProcessors():
                    self.solidSolverTimer.start()
                    self.SolidSolver.run(self.time-self.deltaT, self.time)
                    self.solidSolverTimer.stop()
                    self.solidSolverTimer.cumul()

                # --- Compute and monitor the FSI residual --- #
                res = self.computeSolidInterfaceResidual()
                self.errValue = self.criterion.update(res)
                mpiPrint('\nFSI error value : {}\n'.format(self.errValue), self.mpiComm)
                self.FSIConv = self.criterion.isVerified(self.errValue)

                # --- Initialize d_tilde for the construction of the Wk matrix -- #
                if self.myid in self.manager.getSolidInterfaceProcessors():
                    localSolidInterfaceDisp_X, localSolidInterfaceDisp_Y, localSolidInterfaceDisp_Z = self.SolidSolver.getNodalDisplacements()
                    for iVertex in range(self.manager.getNumberOfLocalSolidInterfaceNodes()):
                        iGlobalVertex = self.manager.getGlobalIndex('solid', self.myid, iVertex)
                        solidInterfaceDisplacement_tilde[iGlobalVertex] = [localSolidInterfaceDisp_X[iVertex], localSolidInterfaceDisp_Y[iVertex], localSolidInterfaceDisp_Z[iVertex]]

                solidInterfaceDisplacement_tilde.assemble()
                
                if self.makeBGS:
                    # --- Relax the solid position --- #
                    mpiPrint('\nProcessing interface displacements...\n', self.mpiComm)
                    self.relaxSolidPosition()
                    self.makeBGS = False
                else:
                    # --- Construct Vk and Wk matrices for the computation of the approximated tangent matrix --- #
                    mpiPrint('\nCorrect solid interface displacements using IQN-MVJ method...\n', self.mpiComm)
                    
                    # --- Start gathering on root process --- #
                    res_X_Gat, res_Y_Gat, res_Z_Gat = mpiGatherInterfaceData(res, ns+d, self.mpiComm, 0)
                    solidInterfaceResidual0_X_Gat, solidInterfaceResidual0_Y_Gat, solidInterfaceResidual0_Z_Gat = mpiGatherInterfaceData(solidInterfaceResidual0, ns+d, self.mpiComm, 0)
                    solidInterfaceDisplacement_tilde_X_Gat, solidInterfaceDisplacement_tilde_Y_Gat, solidInterfaceDisplacement_tilde_Z_Gat = mpiGatherInterfaceData(solidInterfaceDisplacement_tilde, ns, self.mpiComm, 0)
                    solidInterfaceDisplacement_tilde1_X_Gat, solidInterfaceDisplacement_tilde1_Y_Gat, solidInterfaceDisplacement_tilde1_Z_Gat = mpiGatherInterfaceData(solidInterfaceDisplacement_tilde1, ns, self.mpiComm, 0)
                    if self.myid == 0:
                        res_X_Gat_C = res_X_Gat[:ns] # Copies for operating on, length=ns, not ns+d
                        res_Y_Gat_C = res_Y_Gat[:ns]
                        res_Z_Gat_C = res_Z_Gat[:ns]
                        solidInterfaceResidual0_X_Gat_C = solidInterfaceResidual0_X_Gat[:ns] # Copies for operating on, length=ns, not ns+d
                        solidInterfaceResidual0_Y_Gat_C = solidInterfaceResidual0_Y_Gat[:ns]
                        solidInterfaceResidual0_Z_Gat_C = solidInterfaceResidual0_Z_Gat[:ns]

                        if self.FSIIter == 0: # Use J^-1 from the previous time step because Vk and Wk are empty
                            
                            invJ = self.invJprev.copy()

                        else: # Vk and Wk matrices are enriched only starting from the second iteration of every FSI loop
                            if self.manager.nDim == 3:
                                delta_res = np.concatenate([res_X_Gat_C - solidInterfaceResidual0_X_Gat_C, res_Y_Gat_C - solidInterfaceResidual0_Y_Gat_C, res_Z_Gat_C - solidInterfaceResidual0_Z_Gat_C], axis=0)
                                delta_d = np.concatenate([solidInterfaceDisplacement_tilde_X_Gat - solidInterfaceDisplacement_tilde1_X_Gat, solidInterfaceDisplacement_tilde_Y_Gat - solidInterfaceDisplacement_tilde1_Y_Gat, solidInterfaceDisplacement_tilde_Z_Gat - solidInterfaceDisplacement_tilde1_Z_Gat], axis = 0)
                            else:
                                delta_res = np.concatenate([res_X_Gat_C - solidInterfaceResidual0_X_Gat_C, res_Y_Gat_C - solidInterfaceResidual0_Y_Gat_C], axis=0)
                                delta_d = np.concatenate([solidInterfaceDisplacement_tilde_X_Gat - solidInterfaceDisplacement_tilde1_X_Gat, solidInterfaceDisplacement_tilde_Y_Gat - solidInterfaceDisplacement_tilde1_Y_Gat], axis = 0)
                            
                            Vk.insert(0, delta_res)
                            Wk.insert(0, delta_d)
                            
                            nIt+=1
                        
                            Vk_mat = np.vstack(Vk).T
                            Wk_mat = np.vstack(Wk).T
                            
                            Vinv = np.linalg.pinv(Vk_mat)
                            self.invJ = self.invJprev+np.dot(Wk_mat-self.invJprev.dot(Vk_mat),Vinv)
                            invJ = self.invJ.copy()

                        if self.manager.nDim == 3:
                            Res = np.concatenate([res_X_Gat_C, res_Y_Gat_C, res_Z_Gat_C], axis=0)
                        else:
                            Res = np.concatenate([res_X_Gat_C, res_Y_Gat_C], axis=0)

                        np.fill_diagonal(invJ,invJ.diagonal()-1)

                        if self.manager.nDim == 3:
                            delta_ds_loc = np.split(np.dot(invJ,-Res),3,axis=0)
                            
                            delta_ds_loc_X = delta_ds_loc[0]
                            delta_ds_loc_Y = delta_ds_loc[1]
                            delta_ds_loc_Z = delta_ds_loc[2]
                        else:
                            delta_ds_loc = np.split(np.dot(invJ,-Res),2,axis=0)
                            
                            delta_ds_loc_X = delta_ds_loc[0]
                            delta_ds_loc_Y = delta_ds_loc[1]
                            delta_ds_loc_Z = np.zeros(ns)
                        
                        for iVertex in range(delta_ds_loc_X.shape[0]):
                            iGlobalVertex = self.manager.getGlobalIndex('solid', self.myid, iVertex)
                            delta_ds[iGlobalVertex] = [delta_ds_loc_X[iVertex], delta_ds_loc_Y[iVertex], delta_ds_loc_Z[iVertex]]
                    
                    # --- Go back to parallel run --- #
                    mpiBarrier(self.mpiComm)
                    delta_ds.assemble()
                    self.interfaceInterpolator.solidInterfaceDisplacement += delta_ds
                
                if self.computeTangentMatrixBasedOnFirstIt:
                    if self.FSIIter == 0:
                        res.copy(solidInterfaceResidual0)
                        solidInterfaceDisplacement_tilde.copy(solidInterfaceDisplacement_tilde1)
                else:
                    res.copy(solidInterfaceResidual0)
                    solidInterfaceDisplacement_tilde.copy(solidInterfaceDisplacement_tilde1)
            
            if self.writeInFSIloop == True:
                self.writeRealTimeData()
            
            self.FSIIter += 1
            
        # --- Update the FSI history file --- #
        if self.timeIter > self.timeIterTreshold:
            mpiPrint('\n*************** IQN-MVJ is converged ***************', self.mpiComm)