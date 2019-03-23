#!/usr/bin/env python
# -*- coding: latin-1; -*-

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

ModalInterface.py
Python interface between the modal data base and CUPyDO.
Authors 

'''

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import numpy as np
from cupydo.genericSolvers import SolidSolver
from math import *
from scipy import integrate
from numpy.linalg import inv

# ----------------------------------------------------------------------
#  Modal solver interface class
# ----------------------------------------------------------------------

class modalInterpreter(SolidSolver):
    """
    Description
    """

    def __init__(self, confFile, initialModalDisp, initialModalVel, modalMass, modalDamping, modalStiffness):
        """
        Des.
        """

        print("\n***************************** Initialize modal interpreter *****************************\n")

        self.modalDataFileName = confFile       
        self.interfaceID = None
        self.nNodes = self.__linecount(self.modalDataFileName)-1
        self.nHaloNode = 0
        self.nPhysicalNodes = self.nNodes

        self.nodalGlobalIndex = np.zeros(self.nPhysicalNodes, dtype=int)
                
        self.haloNodeList = {}

        self.nodalCoord_X = np.zeros(self.nPhysicalNodes, dtype=float)
        self.nodalCoord_Y = np.zeros(self.nPhysicalNodes, dtype=float)
        self.nodalCoord_Z = np.zeros(self.nPhysicalNodes, dtype=float)
        
        self.y0 = np.concatenate((initialModalDisp, initialModalVel))
        self.Mq = modalMass
        self.invMq = inv(self.Mq)
        self.Kq = modalStiffness
        self.Cq = modalDamping
        
        self.nModes = len(initialModalDisp)
        self.nodalMod_X = np.zeros((self.nPhysicalNodes, self.nModes), dtype=float)
        self.nodalMod_Y = np.zeros((self.nPhysicalNodes, self.nModes), dtype=float)
        self.nodalMod_Z = np.zeros((self.nPhysicalNodes, self.nModes), dtype=float)
        
        self.__readModalData()
               
        self.Phi = np.concatenate((self.nodalMod_X, self.nodalMod_Y, self.nodalMod_Z))            
        self.PhiT = self.Phi.transpose()
        # print(np.amax(self.Phi[:, 0]))
        # print(np.amax(self.Phi[:, 1]))
        # print(np.amax(self.Phi[:, 2]))
        # print(np.amax(self.Phi[:, 3]))
        # print(np.amax(self.Phi[:, 4]))
        # print(np.amax(self.Phi[:, 5]))
        
        self.x = np.dot(self.Phi, self.y0[0:self.nModes])
        
        self.nodalDisp_X = self.x[0:self.nPhysicalNodes]
        self.nodalDisp_Y = self.x[self.nPhysicalNodes:2*self.nPhysicalNodes]
        self.nodalDisp_Z = self.x[2*self.nPhysicalNodes:3*self.nPhysicalNodes]
        
        self.nodalLoad_X = np.zeros(self.nPhysicalNodes, dtype=float)
        self.nodalLoad_Y = np.zeros(self.nPhysicalNodes, dtype=float)
        self.nodalLoad_Z = np.zeros(self.nPhysicalNodes, dtype=float)      

        self.f = np.concatenate((self.nodalLoad_X, self.nodalLoad_Y, self.nodalLoad_Z))
        self.fq = np.dot(self.PhiT, self.f)
        
        f = open('output_t_y0_fq.txt', 'w')
        f.write('Time (s), ')
        for i in range(len(self.y0)):
            f.write('y0_{}, '.format(i + 1))
        for i in range(len(self.fq)):
            if i == len(self.fq) - 1:
                f.write('fq_{}\n'.format(i + 1))
            else:
                f.write('fq_{}, '.format(i + 1))
            
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
                    for iMode in range(self.nModes):
                        self.nodalMod_X[iVertex, iMode] = float(line[4 + 3*iMode])
                        self.nodalMod_Y[iVertex, iMode] = float(line[5 + 3*iMode])
                        self.nodalMod_Z[iVertex, iMode] = float(line[6 + 3*iMode])
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
        
        def f(t, y, self):           
                        
            return np.concatenate([y[self.nModes:2*self.nModes], np.dot(self.invMq, (-np.dot(self.Cq, y[self.nModes:2*self.nModes]) - np.dot(self.Kq, y[0:self.nModes]) + self.fq))])

        if t2 > 0:
            t = np.array([t1, t2])
            y = np.zeros((len(t), len(self.y0)))
            y[0, :] = self.y0
            
            r = integrate.ode(f).set_integrator("dopri5") # Explicit runge-kutta method of order (4)5 due to Dormand & Prince
            r.set_initial_value(self.y0, t1).set_f_params(self)
            for i in range(1, len(t)):
               y[i, :] = r.integrate(t[i])
               if not r.successful():
                   raise RuntimeError("Could not integrate")
                   
            self.y0 = y[1, :]
                          
            self.x = np.dot(self.Phi, self.y0[0:self.nModes])
            
            self.nodalDisp_X = self.x[0:self.nPhysicalNodes]
            self.nodalDisp_Y = self.x[self.nPhysicalNodes:2*self.nPhysicalNodes]
            self.nodalDisp_Z = self.x[2*self.nPhysicalNodes:3*self.nPhysicalNodes]      
           
    def applyNodalLoads(self, load_X, load_Y, load_Z, time):
        
        self.nodalLoad_X = load_X
        self.nodalLoad_Y = load_Y
        self.nodalLoad_Z = load_Z
        
        self.f = np.concatenate((self.nodalLoad_X, self.nodalLoad_Y, self.nodalLoad_Z))
        self.fq = np.dot(self.PhiT, self.f)
        
        f = open('output_t_y0_fq.txt', 'a')
        f.write('{0:.8e}, '.format(time))
        for i in range(len(self.y0)):
            f.write('{0:.8e}, '.format(self.y0[i]))
        for i in range(len(self.fq)):
            if i == len(self.fq) - 1:
                f.write('{0:.8e}\n'.format(self.fq[i]))
            else:
                f.write('{0:.8e}, '.format(self.fq[i]))
    
        return            
            
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


