#! /usr/bin/env python
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

interfaceData.py
Matrix and vector representation of interface data.
Authors : David THOMAS, Marco Lucio CERQUAGLIA, Romain BOMAN

'''

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import numpy as np
import scipy as sp
import sys

import ccupydo

np.set_printoptions(threshold=sys.maxsize)

# ----------------------------------------------------------------------
#   FlexInterfaceData class
# ----------------------------------------------------------------------

class FlexInterfaceData(ccupydo.CFlexInterfaceData):
    """
    Description
    """

    def __init__(self, val_nPoint, val_nDim, mpiComm=None):
        """
        Des.
        """
        self.mpiComm = mpiComm
        
        ccupydo.CFlexInterfaceData.__init__(self, val_nPoint, val_nDim, mpiComm)

        #self.mpiComm = mpiComm
        #self.nPoint = val_nPoint
        #self.nDim = val_nDim

        #self.dataContainer = []

        #if mpiComm != None:
        #    for iDim in range(self.nDim):
        #        from petsc4py import PETSc
        #        data = PETSc.Vec().create(self.mpiComm)
        #        data.setType('mpi')
        #        data.setSizes(self.nPoint)
        #        data.set(0.0)
        #        self.dataContainer.append(data)
        #    self.myid = self.mpiComm.Get_rank()
        #    self.mpiSize = self.mpiComm.Get_size()
        #    startIndex , stopIndex = self.dataContainer[0].getOwnershipRange()
        #    #!!! stopIndex is 1 more than the true index !!!
        #    #startIndex , stopIndex = self.getOwnershipRange()
        #    self.indexing = self.mpiComm.allgather((startIndex , stopIndex))
        #else:
        #    for iDim in range(self.nDim):
        #        data = np.zeros(self.nPoint, dtype=float)
        #        self.dataContainer.append(data)
        #    self.myid = 0
        #    self.mpiSize = 1

    def __setitem__(self, index, values):
        """
        Des.
        """

        if type(values) != list:
            raise TypeError("FlexInterfaceData.__setitem__ needs list as argument !")

        if len(values) != self.nDim:
            raise IndexError("Length of values does not match nDim !")
        else:
            for iDim in range(self.nDim):
                self.setValue(iDim, index, values[iDim])

    def __add__(self, dataToAdd):
        """
        Des.
        """

        if type(dataToAdd) == type(self):
            if self.nDim != dataToAdd.nDim:
                raise IndexError("Dimensions do not match for + operator !")
            if self.nPoint != dataToAdd.nPoint:
                raise IndexError("Lengthes do not match for + operator !")

        newData = FlexInterfaceData(self.nPoint, self.nDim, self.comm)
        self.copy(newData)
        newData.add(dataToAdd)

        return newData

    def __radd__(self, dataToAdd):
        """
        Des.
        """

        newData = self + dataToAdd

        return newData

    def __iadd__(self, dataToAdd):
        """
        Des.
        """

        if type(dataToAdd) == type(self):
            if self.nDim != dataToAdd.nDim:
                raise IndexError("Dimensions do not match for + operator !")
            if self.nPoint != dataToAdd.nPoint:
                raise IndexError("Lengthes do not match for + operator !")

        self.add(dataToAdd)

        return self

    def __sub__(self, dataToSub):
        """
        Des.
        """

        if type(dataToSub) == type(self):
            if self.nDim != dataToSub.nDim:
                raise IndexError("Dimensions do not match for + operator !")
            if self.nPoint != dataToSub.nPoint:
                raise IndexError("Lengthes do not match for + operator !")

        newData = FlexInterfaceData(self.nPoint, self.nDim, self.comm)
        self.copy(newData)
        newData.sub(dataToSub)

        return newData

    def __rsub__(self, dataToSub):
        """
        Des.
        """

        if type(dataToSub) == type(self):
            if self.nDim != dataToSub.nDim:
                raise IndexError("Dimensions do not match for + operator !")
            if self.nPoint != dataToSub.nPoint:
                raise IndexError("Lengthes do not match for + operator !")

        newData = -1*self + dataToSub

        return newData

    def __isub__(self, dataToSub):
        """
        Des.
        """

        if type(dataToSub) == type(self):
            if self.nDim != dataToSub.nDim:
                raise IndexError("Dimensions do not match for + operator !")
            if self.nPoint != dataToSub.nPoint:
                raise IndexError("Lengthes do not match for + operator !")

        self.sub(dataToSub)

        return self

    def __mul__(self, mulVal):
        """
        Des.
        """

        newData = FlexInterfaceData(self.nPoint, self.nDim, self.comm)
        self.copy(newData)
        newData.scale(mulVal)

        return newData

    def __rmul__(self, mulVal):
        """
        Des
        """

        newData = self*mulVal

        return newData

    def __imul__(self, mulVal):
        """
        Des.
        """

        self.scale(mulVal)

        return self
    
    def dot(self, dataToDot):
        
        dotList = []

        if self.mpiComm != None:
            dotList = ccupydo.CFlexInterfaceData.dot(self, dataToDot)
        else:
            for iDim in range(self.nDim):
                myData = self.getData(iDim)
                dotData = dataToDot.getData(iDim)
                val_dot = myData.dot(dotData)
                dotList.append(val_dot)

        return dotList
    
    def sum(self):
        
        sumList = []

        if self.mpiComm != None:
            sumList = ccupydo.CFlexInterfaceData.sum(self)
        else:
            for iDim in range(self.nDim):
                myData = self.getData(iDim)
                val_sum = myData.sum()
                sumList.append(val_sum)

        return sumList
    
    def norm(self):
        
        normList = []

        if self.mpiComm != None:
            normList = ccupydo.CFlexInterfaceData.norm(self)
        else:
            for iDim in range(self.nDim):
                myData = self.getData(iDim)
                val_norm = np.linalg.norm(myData)
                normList.append(val_norm)

        return normList

# ----------------------------------------------------------------------
#    InterfaceMatrix class
# ----------------------------------------------------------------------

class InterfaceMatrix(ccupydo.CInterfaceMatrix):
    """
    Define a matrix based on fluid-structure interface data.
    Designed for parallel computations (also works in serial).
    Inherited public members :
        -createDense()
        -createSparse()
        -createSparseFullAlloc()
        -setValue()
        -assemble()
        -getMat()
    """

    def __init__(self, sizes, mpiComm=None):
        """
        Overloaded constructor
        """

        ccupydo.CInterfaceMatrix.__init__(self, sizes[0],sizes[1])

        self.sizes = sizes
        self.mpiComm = mpiComm

    def mult(self, Data , DataOut):
        """
        Performs interface matrix-data multiplication.
        """

        if self.mpiComm != None:
            ccupydo.CInterfaceMatrix.mult(self, Data, DataOut)
        else:
            PyH = self.getMat();
            dim = Data.getDim()
            for iDim in range(dim):
                np.dot(PyH, Data.getData(iDim), DataOut.getData(iDim))
