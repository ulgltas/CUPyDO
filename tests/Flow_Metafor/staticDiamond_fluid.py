#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# CUPyDO configuration file for Flow
# Diamond airfoil
# Adrien Crovato

import os

def getParams():
    p = {}
    # Input/Output
    p['File'] = os.path.abspath(os.path.join(os.path.dirname(__file__), 'models/diamond_fluid.geo')) # Input file containing the model
    p['Pars'] = {'xLgt' : 5, 'yLgt' : 5, 'msF' : 1., 'msA' : 0.01} # Parameters for input file model
    p['Dim'] = 2 # Problem dimension
    p['Format'] = 'gmsh' # Save format (vtk or gmsh)
    # Markers
    p['Fluid'] = 'field' # Name of physical group containing the fluid
    p['Farfield'] = ['upstream', 'side', 'downstream'] # LIST of name of physical groups containing the farfield boundaries (upstream/downstream should be first/last element)
    p['Wing'] = 'airfoil' # Name of physical group containing the airfoil FSI boundary
    p['Wake'] = 'wake' # Name of physical group containing the wake
    p['Te'] =  'te' # Name of physical group containing the trailing edge
    # Freestream
    p['M_inf'] = 0.31 # Freestream Mach number
    p['AoA'] = 3. # Freestream angle of attack
    p['P_dyn'] = 0.5*102*102 # Dynamic pressure
    # Geometry
    p['S_ref'] = 1. # Reference surface length (c_ref for 2D)
    p['c_ref'] = 1. # Reference chord length
    p['x_ref'] = 0. # Reference point for moment computation (x)
    p['y_ref'] = 0. # Reference point for moment computation (y)
    p['z_ref'] = 0. # Reference point for moment computation (z)
    # Numerical
    p['LSolver'] = 'Pardiso' # Linear solver (Pardiso, GMRES, MUMPS or SparseLU)
    p['NSolver'] = 'Newton' # Noninear solver type (Picard or Newton)
    p['Rel_tol'] = 1e-6 # Relative tolerance on solver residual
    p['Abs_tol'] = 1e-8 # Absolute tolerance on solver residual
    p['Max_it'] = 25 # Solver maximum number of iterations
    p['LS_tol'] = 1e-6 # Tolerance on linesearch residual
    p['Max_it_LS'] = 10 # Linesearch maximum number of iterations
    p['AV_thrsh'] = 1e-2 # Residual threshold below which the artificial viscosity is decreased
    p['M_crit'] = 5. # Critical Mach number above which density is damped
    return p
