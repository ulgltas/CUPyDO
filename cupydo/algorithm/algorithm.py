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
from ..timeStep import TimeStep

np.set_printoptions(threshold=sys.maxsize)

# ----------------------------------------------------------------------
#    Algorithm class
# ----------------------------------------------------------------------

class Algorithm(object):
    def __init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, p, mpiComm):

        mpiPrint("Initializing FSI Algorithm",mpiComm,titlePrint)

        self.mpiComm = mpiComm
        self.manager = Manager
        self.FluidSolver = FluidSolver
        self.SolidSolver = SolidSolver
        self.interfaceInterpolator = InterfaceInterpolator
        
        self.globalTimer = Timer()
        self.meshDefTimer = Timer()
        self.communicationTimer = Timer()
        self.fluidSolverTimer = Timer()
        self.solidSolverTimer = Timer()
        self.solidUpdateTimer = Timer()
        self.fluidUpdateTimer = Timer()

        self.totTime = p['tTot']
        self.interpType = p['interpType']
        self.step = TimeStep(Manager, FluidSolver, SolidSolver, p, mpiComm)
        
        if self.mpiComm != None:
            self.myid = self.mpiComm.Get_rank()
            self.mpiSize = self.mpiComm.Get_size()
        else:
            self.myid = 0
            self.mpiSize = 1

        self.solidHasRun = False
        self.verified = True

    # --- Function called when the FSI coupling failed --- #
    def resetInternalVars(self): return

    def setFSIInitialConditions(self):

        if self.manager.mechanical:
            self.interfaceInterpolator.getDisplacementFromSolidSolver()
            self.solidToFluidMechaTransfer()
            self.FluidSolver.setInitialMeshDeformation()
        if self.manager.thermal:
            self.interfaceInterpolator.getTemperatureFromSolidSolver()
            self.solidToFluidThermalTransfer()
            
        self.FluidSolver.boundaryConditionsUpdate()

    def solidToFluidMechaTransfer(self):

        self.communicationTimer.start()
        self.interfaceInterpolator.interpolateSolidDisplacementOnFluidMesh()
        self.interfaceInterpolator.setDisplacementToFluidSolver(self.step.dt)
        self.communicationTimer.stop()
        self.communicationTimer.cumul()

    def solidToFluidThermalTransfer(self):

        self.communicationTimer.start()
        if self.interfaceInterpolator.chtTransferMethod == 'TFFB' or self.interfaceInterpolator.chtTransferMethod == 'hFFB':
            self.interfaceInterpolator.interpolateSolidHeatFluxOnFluidMesh()
            self.interfaceInterpolator.setHeatFluxToFluidSolver(self.step.dt)
        elif self.interfaceInterpolator.chtTransferMethod == 'FFTB' or self.interfaceInterpolator.chtTransferMethod == 'hFTB':
            self.interfaceInterpolator.interpolateSolidTemperatureOnFluidMesh()
            self.interfaceInterpolator.setTemperatureToFluidSolver(self.step.dt)
        else: raise Exception('Wrong CHT transfer method, use: TFFB, FFTB, hFTB, hFFB')
        self.communicationTimer.stop()
        self.communicationTimer.cumul()

    def fluidToSolidMechaTransfer(self):

        self.communicationTimer.start()
        if self.interpType == 'conservative':
            self.interfaceInterpolator.getForceFromFluidSolver()
            self.interfaceInterpolator.interpolateFluidLoadsOnSolidMesh()
            self.interfaceInterpolator.setForceToSolidSolver(self.step.dt)
        elif self.interpType == 'consistent':
            self.interfaceInterpolator.getStressFromFluidSolver()
            self.interfaceInterpolator.interpolateFluidLoadsOnSolidMesh()
            self.interfaceInterpolator.setStressToSolidSolver(self.step.dt)
        else: raise Exception('Wrong interpolation type, use: conservative, consistent')
        self.communicationTimer.stop()
        self.communicationTimer.cumul()
        
    def fluidToSolidThermalTransfer(self):

        self.communicationTimer.start()
        if self.interfaceInterpolator.chtTransferMethod == 'TFFB':
            self.interfaceInterpolator.getTemperatureFromFluidSolver()
            self.interfaceInterpolator.interpolateFluidTemperatureOnSolidMesh()
            self.interfaceInterpolator.setTemperatureToSolidSolver(self.step.dt)
        elif self.interfaceInterpolator.chtTransferMethod == 'FFTB':
            self.interfaceInterpolator.getHeatFluxFromFluidSolver()
            self.interfaceInterpolator.interpolateFluidHeatFluxOnSolidMesh()
            self.interfaceInterpolator.setHeatFluxToSolidSolver(self.step.dt)
        elif self.interfaceInterpolator.chtTransferMethod == 'hFFB' or self.interfaceInterpolator.chtTransferMethod == 'hFTB':
            self.interfaceInterpolator.getRobinTemperatureFromFluidSolver()
            self.interfaceInterpolator.interpolateFluidRobinTemperatureOnSolidMesh()
            self.interfaceInterpolator.setRobinHeatFluxToSolidSolver(self.step.dt)
        else: raise Exception('Wrong CHT transfer method, use: TFFB, FFTB, hFTB, hFFB')
        self.communicationTimer.stop()
        self.communicationTimer.cumul()

    def fluidToSolidAdjointTransfer(self):

        self.communicationTimer.start()
        self.interfaceInterpolator.getAdjointDisplacementFromFluidSolver()
        self.interfaceInterpolator.interpolateFluidAdjointDisplacementOnSolidMesh()
        self.interfaceInterpolator.setAdjointDisplacementToSolidSolver(self.step.dt)
        self.communicationTimer.stop()
        self.communicationTimer.cumul()

    def solidToFluidAdjointTransfer(self):

        self.communicationTimer.start()
        if self.interpType == 'conservative':
            self.interfaceInterpolator.getAdjointForceFromSolidSolver()
            self.interfaceInterpolator.interpolateSolidAdjointLoadsOnFluidMesh()
            self.interfaceInterpolator.setAdjointForceToFluidSolver(self.step.dt)
        elif self.interpType == 'consistent':
            self.interfaceInterpolator.getAdjointStressFromSolidSolver()
            self.interfaceInterpolator.interpolateSolidAdjointLoadsOnFluidMesh()
            self.interfaceInterpolator.setAdjointStressToFluidSolver(self.step.dt)
        self.communicationTimer.stop()
        self.communicationTimer.cumul()

# ----------------------------------------------------------------------
#    Explicit Algorithm class
# ----------------------------------------------------------------------

class AlgorithmExplicit(Algorithm):
    def __init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, p, mpiComm):
        Algorithm.__init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, p, mpiComm)

    def run(self):
        
        # --- Initialize the algorithm --- #
        mpiPrint("Begin FSI Computation",self.mpiComm,titlePrint)
        self.initInterfaceData()
        self.iniRealTimeData()
        self.globalTimer.start()
        self.setFSIInitialConditions()
        
        try:
            if self.manager.regime == 'unsteady':
                self.__unsteadyRun()
            else:
                raise Exception('Explicit coupling is only valid for unsteady computations!')
        except:
            mpiPrint('\nError when executing Explicit coupling\n', self.mpiComm)
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

            # --- Preprocess the temporal iteration --- #
            self.FluidSolver.preprocessTimeIter(self.step.timeIter)
            if self.myid in self.manager.getSolidSolverProcessors():
                self.SolidSolver.preprocessTimeIter(self.step.timeIter)

            # --- Internal FSI loop --- #
            self.criterion.reset()
            self.verified = self.fsiCoupling()
            self.totNbOfFSIIt += self.FSIIter
            mpiBarrier(self.mpiComm)

            # --- Update TimeStep class and restart if FSI failed --- #
            if not self.verified:
                mpiPrint('\nFSI coupling did not converge, restart with smaller time step\n', self.mpiComm)
                self.step.updateTime(self.verified)
                self.resetInternalVars()
                self.writeRealTimeData()
                continue

            # --- Update the fluid and solid solver for the next time step --- #
            if self.myid in self.manager.getSolidSolverProcessors():
                self.solidUpdateTimer.start()
                self.SolidSolver.update()
                self.solidUpdateTimer.stop()
                self.solidUpdateTimer.cumul()

            self.fluidUpdateTimer.start()
            self.FluidSolver.update(self.step.dt)
            self.fluidUpdateTimer.stop()
            self.fluidUpdateTimer.cumul()

            # --- Update TimeStep class, export the results and write FSI history --- #
            self.step.updateTime(self.verified)
            self.step.updateSave()
            self.writeRealTimeData()

        # --- End of the temporal loop --- #

    def fsiCoupling(self):
        """
        FSI explicit coupling
        """

        mpiPrint("Enter Explicit FSI Coupling",self.mpiComm,titlePrint)

        # --- Solid to fluid mechanical transfer --- #
        if self.manager.mechanical:
            self.solidToFluidMechaTransfer()
            mpiPrint('\nPerforming mesh deformation...\n', self.mpiComm)
            self.meshDefTimer.start()
            self.FluidSolver.meshUpdate(self.step.timeIter)
            self.meshDefTimer.stop()
            self.meshDefTimer.cumul()

        # --- Solid to fluid thermal transfer --- #
        if self.manager.thermal:
            self.solidToFluidThermalTransfer()

        self.FluidSolver.boundaryConditionsUpdate()

        # --- Fluid solver call  --- #
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
        self.solidHasRun = True

        # --- Check if the solid solver succeeded --- #
        try: solidProc = int(self.manager.getSolidSolverProcessors())
        except: raise Exception('Only one solid solver process is supported yet')
        verif = mpiScatter(verif, self.mpiComm, solidProc)
        self.solidHasRun = True
        if not verif: return False

        return True

    def iniRealTimeData(self):

        if self.myid in self.manager.getSolidSolverProcessors():
            self.SolidSolver.initRealTimeData()

    def writeRealTimeData(self):

        if self.myid == 0:
            self.FluidSolver.saveRealTimeData(self.step.time, 0)
            self.SolidSolver.saveRealTimeData(self.step.time, 0)

    def printExitInfo(self):
        
        mpiPrint('[cpu FSI total]: ' + str(self.globalTimer.cumulTime) + ' s', self.mpiComm)
        mpiPrint('[cpu FSI fluid mesh mapping]: ' + str(self.interfaceInterpolator.mappingTimer.cumulTime) + ' s', self.mpiComm)
        mpiPrint('[cpu FSI fluid mesh deformation]: ' + str(self.meshDefTimer.cumulTime) + ' s', self.mpiComm)
        mpiPrint('[cpu FSI communications]: ' + str(self.communicationTimer.cumulTime) + ' s', self.mpiComm)
        mpiPrint('[cpu FSI fluid solver]: ' + str(self.fluidSolverTimer.cumulTime) + ' s', self.mpiComm)
        mpiPrint('[cpu FSI solid solver]: ' + str(self.solidSolverTimer.cumulTime) + ' s', self.mpiComm)
        mpiPrint('[cpu FSI fluid remeshing]: ' + str(self.fluidUpdateTimer.cumulTime) + ' s', self.mpiComm)
        mpiPrint('[cpu FSI solid remeshing]: ' + str(self.solidUpdateTimer.cumulTime) + ' s', self.mpiComm)
        mpiPrint('[Time steps FSI]: ' + str(self.step.timeIter), self.mpiComm)
        mpiPrint('[Successful Run FSI]: ' + str(self.step.time >= (self.totTime - 2*self.step.dt)), self.mpiComm) # NB: self.totTime - 2*self.step.dt is the extreme case that can be encountered due to rounding effects!
        mpiPrint('[Mean n. of FSI Iterations]: ' + str(1), self.mpiComm)

        if self.myid == 0 :
            self.FluidSolver.printRealTimeData(self.step.time, 0)
            self.SolidSolver.printRealTimeData(self.step.time, 0)

        mpiPrint('RES-FSI-FSIhistory: ' + str(self.step.timeIter) + '\t' + str(self.step.time) + '\t' + str(1.0) + '\t' + str(1) + '\n', self.mpiComm)
