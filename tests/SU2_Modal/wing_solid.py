#! /usr/bin/env python
# -*- coding: utf-8; -*-

import numpy as np
import ModalSolver as ms

class Module:
    def __init__(self, _solver):
        self.solver = _solver

def getModal():
    # Config file
    fconfig = '../../tests/SU2_Modal/models/wing_solid.csv'
    # Paramters
    initialModalDisp = np.array([5.6023, 0., 0., 0., 0., 0.])
    initialModalVel = np.zeros(6, dtype=float)
    modalMass = np.diag([1., 1., 1., 1., 1., 1.])
    modalDamping = np.zeros((6, 6), dtype=float)
    modalStiffness = np.diag([1.817819061E+02, 1.178750000E+03, 1.709068970E+03, 5.454931152E+03, 8.660587891E+03, 9.815511719E+03])
    nModes = initialModalDisp.shape[0]

    # Solver
    solver = ms.ModalSolver()
    solver.setMatrices(nModes, modalMass, modalDamping, modalStiffness)
    solver.readModes(fconfig)
    solver.setInitial(initialModalDisp, initialModalVel)

    return Module(solver)