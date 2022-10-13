#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# CUPyDO configuration file for DART
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
    p['Wings'] = ['wing'] # LIST of names of physical group containing the lifting surface boundary (first element will be the FSI boundary)
    p['Wakes'] = ['wake'] # LIST of names of physical groups containing the wake
    p['WakeTips'] = ['wakeTip'] # LIST of names of physical groups containing the edge of the wake
    p['Tes'] = ['te'] # LIST of names of physical groups containing the trailing edge
    # Freestream
    p['M_inf'] = 0.8 # Freestream Mach number
    p['AoA'] = 1. # Freestream angle of attack
    p['Q_inf'] = 0.5*0.094*247.1*247.1 # Dynamic pressure
    # Geometry
    p['S_ref'] = .35 # Reference surface length (c_ref for 2D)
    p['c_ref'] = .47 # Reference chord length
    p['x_ref'] = 0. # Reference point for moment computation (x)
    p['y_ref'] = 0. # Reference point for moment computation (y)
    p['z_ref'] = 0. # Reference point for moment computation (z)
    # Numerical
    p['LSolver'] = 'GMRES' # Linear solver (PARDISO, GMRES, or MUMPS)
    p['G_fill'] = 2 # Fill-in factor for GMRES preconditioner
    p['G_tol'] = 1e-5 # Tolerance for GMRES
    p['G_restart'] = 50 # Restart for GMRES
    p['Rel_tol'] = 1e-6 # Relative tolerance on solver residual
    p['Abs_tol'] = 1e-8 # Absolute tolerance on solver residual
    p['Max_it'] = 25 # Solver maximum number of iterations
    return p
