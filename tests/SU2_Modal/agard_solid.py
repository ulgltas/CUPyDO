#! /usr/bin/env python
# -*- coding: utf-8; -*-

import numpy as np
import ModalSolver as ms
import os.path

class Module:
    def __init__(self, _solver):
        self.solver = _solver

def getModal():
    # Config file
    fconfig = os.path.join(os.path.abspath(os.path.dirname(__file__)),'models/agard_modes.csv')
    # Paramters
    initialModalDisp = np.zeros(4, dtype=float)
    initialModalVel = np.zeros(4, dtype=float)
    initialModalLoads = np.zeros(4, dtype=float)
    modalMass = np.diag([2.9107e-04, 8.3181e-05, 1.7447e-04, 3.4281e-05])
    modalDamping = np.zeros((4, 4), dtype=float)
    modalStiffness = np.diag([1.0468, 5.3468, 17.3717, 12.9114])
    nModes = initialModalDisp.shape[0]

    # Solver
    solver = ms.ModalSolver(nModes)
    solver.setMatrices(modalMass, modalDamping, modalStiffness)
    solver.readModes(fconfig)
    solver.setInitial(initialModalDisp, initialModalVel, initialModalLoads)
    solver.setExtractor([16, 13808])

    return Module(solver)