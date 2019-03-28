#! /usr/bin/env python
# -*- coding: utf-8; -*-

import numpy as np
import ModalSolver as ms

class Module:
    def __init__(self, _solver):
        self.solver = _solver

def getModal():
    # Config file
    fconfig = 'models/wing_solid.csv'
    # Paramters
    initialModalDisp = np.array([5.6023, 0., 0., 0., 0., 0.])
    initialModalVel = np.zeros(6, dtype=float)
    initialModalForcesX = np.zeros(6, dtype=float)
    initialModalForcesY = np.zeros(6, dtype=float)
    initialModalForcesZ = np.zeros(6, dtype=float)
    modalMass = np.diag([1., 1., 1., 1., 1., 1.])
    modalDamping = np.zeros((6, 6), dtype=float)
    modalStiffness = np.diag([1.817819061E+02, 1.178750000E+03, 1.709068970E+03, 5.454931152E+03, 8.660587891E+03, 9.815511719E+03])
    nModes = initialModalDisp.shape[0]

    # Solver
    solver = ms.ModalSolver()
    solver.setMatrices(nModes, modalMass, modalDamping, modalStiffness)
    solver.setInitial(initialModalDisp, initialModalVel)
    solver.readModes(fconfig)
    solver.updateLoads(initialModalForcesX, initialModalForcesY, initialModalForcesZ)

    return Module(solver)