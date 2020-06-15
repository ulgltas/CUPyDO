#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# CUPyDO configuration file for Flow
# Agard445 wing
# Adrien Crovato

import os

def getParams():
    p = {}
    # Mesh length
    rlems = 0.004 # root leading edge mesh size
    rtems = 0.006 # root trailing edge mesh size
    tlems = 0.002 # tip leading mesh size
    ttems = 0.004 # tip trailing mesh size
    fms = 1.0 # farfield mesh size
    # Input/Output
    p['File'] = os.path.abspath(os.path.join(os.path.dirname(__file__), 'models/agard445_fluid.geo')) # Input file containing the model
    p['Pars'] = {'xL' : 7., 'yL' : 3., 'zL' : 6., 'xO' : -3., 'zO' : -3., 'msLeRt' : rlems, 'msTeRt' : rtems, 'msLeTp' : tlems, 'msTeTp' : ttems, 'msF' : fms} # Parameters for input file model
    p['Dim'] = 3 # Problem dimension
    p['Format'] = 'gmsh' # Save format (vtk or gmsh)
    # Markers
    p['Fluid'] = 'field' # Name of physical group containing the fluid
    p['Symmetry'] = 'symmetry' # Name of physical group containing the (slip) symmetry boundary
    p['Farfield'] = ['upstream', 'farfield', 'downstream'] # LIST of names of physical groups containing the farfield boundaries (upstream/downstream should be first/last element)
    p['Wings'] = ['wing'] # LIST of names of physical group containing the lifting surface boundary
    p['Fsi'] = 'wing' # Name of the physical group containing the FSI boundary
    p['Wakes'] = ['wake'] # LIST of names of physical groups containing the wake
    p['WakeTips'] = ['wakeTip'] # LIST of names of physical groups containing the edge of the wake
    p['TeTips'] = ['teTip'] # LIST of names of physical groups containing the edge of the wake and the trailing edge
    # Freestream
    p['M_inf'] = 0.8 # Freestream Mach number
    p['AoA'] = 1. # Freestream angle of attack
    p['P_dyn'] = 0.5*0.094*247.1*247.1 # Dynamic pressure
    # Geometry
    p['S_ref'] = .35 # Reference surface length (c_ref for 2D)
    p['c_ref'] = .47 # Reference chord length
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
