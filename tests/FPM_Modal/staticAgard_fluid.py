#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# CUPyDO configuration file for Fpm
# Agard445 wing
# Axel Dechamps

import os

def getParams():
    p = {}
    # Mesh length
    nC = 100 # Number of panels along the chord
    nS = 10 # Number of panels along the span
    bC = 1 # Progression along the chord
    pS = 1 # Progression along the span
    # Input/Output
    p['File'] = os.path.abspath(os.path.join(os.path.dirname(__file__), 'models/agard445_fluid.geo')) # Input file containing the model
    p['Pars'] = {'nC' : nC, 'nS' : nS, 'bC' : bC, 'pS' : pS} # Parameters for input file model
    p['Format'] = 'gmsh' # Save format (vtk or gmsh)
    # Markers
    p['Wings'] = ['wing'] # LIST of names of physical group containing the lifting surface boundary
    p['Tes'] = ['te'] # LIST of names of physical groups containing the trailing edges of the lifting surfaces
    p['Fsi'] = 'wing' # Name of the physical group containing the FSI boundary
    # Freestream
    p['M_inf'] = 0. # Freestream Mach number
    p['AoA'] = 1. # Freestream angle of attack [deg]
    p['AoS'] = 0. # Slideslip angle [deg]
    p['P_dyn'] = 0.5*0.094*247.1*247.1 # Dynamic pressure
    # Geometry
    p['LengthWake'] = [10] # LIST of wake lengths corresponding to all the lifting surface boundaries 
    p['S_ref'] = .35 # Reference surface length (c_ref for 2D)
    p['c_ref'] = .47 # Reference chord length
    p['x_ref'] = 0. # Reference point for moment computation (x)
    p['y_ref'] = 0. # Reference point for moment computation (y)
    p['z_ref'] = 0. # Reference point for moment computation (z)
    p['sym'] = 1 # True if computations are done with half a wing
    # Misc
    p['saveFreq'] = 1 # Saving frequency
    return p