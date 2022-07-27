#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# CUPyDO configuration file for ModalSolver
# Agard445 wing
# Adrien Crovato

import os
import numpy as np

def getParams():
    p = {}
    # Input file contining the modes
    p['File'] = os.path.abspath(os.path.join(os.path.dirname(__file__), 'models/agard_modes.csv'))
    # Modal matrices
    p['M_q'] = np.diag([2.9107e-04, 8.3181e-05, 1.7447e-04, 3.4281e-05]) # mass matrix (given as numpy array)
    p['C_q'] = np.zeros((4, 4), dtype=float) # damping matrix (given as numpy array)
    p['K_q'] = np.diag([1.0468, 5.3468, 17.3717, 12.9114]) # stiffness matrix (given as numpy array)
    # Initial conditions
    p['x_i'] = np.zeros(4, dtype=float) # initial modal displacements (given as numpy array)
    p['v_i'] = np.zeros(4, dtype=float) # initial modal velocities (given as numpy array)
    p['f_i'] = np.zeros(4, dtype=float) # initial modal loads (given as numpy array)
    # Number of modes
    p['nm'] = p['x_i'].shape[0]
    # Extractors
    p['Extractors'] = [16, 13808] # python list, cannot be empty
    return p
