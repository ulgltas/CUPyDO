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
    """
    Des.
    """

    def __init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, deltaT, totTime, dtSave, mpiComm=None):
        """
        Des.
        """

        mpiPrint("Initializing FSI Algorithm",mpiComm,titlePrint)

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

        self.totTime = totTime
        self.step = TimeStep(deltaT,dtSave)
        
        if self.mpiComm != None:
            self.myid = self.mpiComm.Get_rank()
            self.mpiSize = self.mpiComm.Get_size()
        else:
            self.myid = 0
            self.mpiSize = 1

        self.solidHasRun = False
        self.verified = True

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
                self.interfaceInterpolator.setDisplacementToFluidSolver(self.step.dt)
                self.FluidSolver.setInitialMeshDeformation()
            else:
                if self.myid in self.manager.getSolidSolverProcessors():
                    self.SolidSolver.setInitialDisplacements()
                self.interfaceInterpolator.getDisplacementFromSolidSolver()
                self.interfaceInterpolator.interpolateSolidDisplacementOnFluidMesh()
                self.interfaceInterpolator.setDisplacementToFluidSolver(self.step.dt)
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
        self.interfaceInterpolator.setLoadsToSolidSolver(self.step.dt)
        self.communicationTimer.stop()
        self.communicationTimer.cumul()

    def solidToFluidMechaTransfer(self):
        """
        Des.
        """

        self.communicationTimer.start()
        self.interfaceInterpolator.interpolateSolidDisplacementOnFluidMesh()
        self.interfaceInterpolator.setDisplacementToFluidSolver(self.step.dt)
        self.communicationTimer.stop()
        self.communicationTimer.cumul()

    def solidToFluidThermalTransfer(self):
        """
        Des.
        """

        self.communicationTimer.start()
        if self.interfaceInterpolator.chtTransferMethod == 'TFFB' or self.interfaceInterpolator.chtTransferMethod == 'hFFB':
            self.interfaceInterpolator.interpolateSolidHeatFluxOnFluidMesh()
            self.interfaceInterpolator.setHeatFluxToFluidSolver(self.step.dt)
        elif self.interfaceInterpolator.chtTransferMethod == 'FFTB' or self.interfaceInterpolator.chtTransferMethod == 'hFTB':
            self.interfaceInterpolator.interpolateSolidTemperatureOnFluidMesh()
            self.interfaceInterpolator.setTemperatureToFluidSolver(self.step.dt)
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
            self.interfaceInterpolator.setTemperatureToSolidSolver(self.step.dt)
        elif self.interfaceInterpolator.chtTransferMethod == 'FFTB':
            self.interfaceInterpolator.getHeatFluxFromFluidSolver()
            self.interfaceInterpolator.interpolateFluidHeatFluxOnSolidMesh()
            self.interfaceInterpolator.setHeatFluxToSolidSolver(self.step.dt)
        elif self.interfaceInterpolator.chtTransferMethod == 'hFFB' or self.interfaceInterpolator.chtTransferMethod == 'hFTB':
            self.interfaceInterpolator.getRobinTemperatureFromFluidSolver()
            self.interfaceInterpolator.interpolateFluidRobinTemperatureOnSolidMesh()
            self.interfaceInterpolator.setRobinHeatFluxToSolidSolver(self.step.dt)
        self.communicationTimer.stop()
        self.communicationTimer.cumul()

# ----------------------------------------------------------------------
#    Explicit Algorithm class
# ----------------------------------------------------------------------

class AlgorithmExplicit(Algorithm):
    """
    Des.
    """

    def __init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, deltaT, totTime, dtSave, mpiComm=None):
        """
        Des.
        """

        Algorithm.__init__(self, Manager, FluidSolver, SolidSolver, InterfaceInterpolator, deltaT, totTime, dtSave, mpiComm)

    def run(self):
        """
        Des.
        """
        
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
                raise Exception('Explicit coupling is only valid for unsteady computations!')
        except:
            mpiPrint('\nA DIVINE ERROR OCCURED...EXITING COMPUTATION\n', self.mpiComm)
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
        """
        Des.
        """

        mpiPrint('Begin time integration\n', self.mpiComm)

        # --- External temporal loop --- #
        while self.step.time < self.totTime:
            
            mpiPrint("\n>>>> Time iteration {} <<<<".format(self.step.timeIter), self.mpiComm)

            # --- Preprocess the temporal iteration --- #
            self.FluidSolver.preprocessTimeIter(self.step.timeIter)
            if self.myid in self.manager.getSolidSolverProcessors():
                self.SolidSolver.preprocessTimeIter(self.step.timeIter)

            # --- Internal FSI loop --- #
            self.fsiCoupling()
            # --- End of FSI loop --- #

            mpiBarrier(self.mpiComm)

            self.totNbOfFSIIt += self.FSIIter

            # --- Update the fluid and solid solver for the next time step --- #
            if self.myid in self.manager.getSolidSolverProcessors():
                self.SolidSolver.update()
            self.FluidSolver.update(self.step.dt)

            # --- Save the fluid and solid solutions ---#


            # --- Update the next required extractor time  ---#


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

            # --- Update the time iteration and FSI history --- #
            self.writeRealTimeData()

        # --- End of the temporal loop --- #

    def fsiCoupling(self):
        """
        FSI explicit coupling
        """

        mpiPrint("Enter Explicit FSI Coupling",self.mpiComm,titlePrint)

        if self.manager.mechanical:
            # --- Solid to fluid mechanical transfer --- #
            self.interfaceInterpolator.getDisplacementFromSolidSolver()
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

        # --- Fluid solver call  --- #
        mpiPrint('\nLaunching fluid solver...', self.mpiComm)
        self.fluidSolverTimer.start()
        self.FluidSolver.run(*self.step.timeFrame())
        self.fluidSolverTimer.stop()
        self.fluidSolverTimer.cumul()
        mpiBarrier(self.mpiComm)

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
            self.SolidSolver.run(*self.step.timeFrame())
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
            self.FluidSolver.saveRealTimeData(self.step.time, 0)
            self.SolidSolver.saveRealTimeData(self.step.time, 0)

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
        mpiPrint('[Time steps FSI]: ' + str(self.step.timeIter), self.mpiComm)
        mpiPrint('[Successful Run FSI]: ' + str(self.step.time >= (self.totTime - 2*self.step.dt)), self.mpiComm) # NB: self.totTime - 2*self.step.dt is the extreme case that can be encountered due to rounding effects!
        mpiPrint('[Mean n. of FSI Iterations]: ' + str(1), self.mpiComm)

        if self.myid == 0 :
            self.FluidSolver.printRealTimeData(self.step.time, 0)
            self.SolidSolver.printRealTimeData(self.step.time, 0)

        mpiPrint('RES-FSI-FSIhistory: ' + str(self.step.timeIter) + '\t' + str(self.step.time) + '\t' + str(1.0) + '\t' + str(1) + '\n', self.mpiComm)