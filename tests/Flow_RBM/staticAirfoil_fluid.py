#! /usr/bin/env python
# -*- coding: utf-8; -*-

"""
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
"""

import math
import flow as f
import tbox
import tbox.gmsh as gmsh

class Module:
    def __init__(self, _msh, _mshDef, _boundary, _solver, _fCp, _dynP):
        # objects
        self.msh = _msh
        self.mshDef = _mshDef
        self.boundary = _boundary
        self.solver = _solver
        self.fCp = _fCp
        # parameters
        self.dynP = _dynP

def initSolver(_solver):
    _solver.relTol = 1e-6
    _solver.absTol = 1e-8
    _solver.lsTol = 1e-6
    _solver.maxIt = 50
    _solver.maxLsIt = 10
    _solver.avThrsh = 1e-3

def getFlow():

    # define flow variables
    alpha = 3*math.pi/180
    U_inf = [math.cos(alpha), math.sin(alpha)] # norm should be 1
    M_inf = 0.03
    gamma = 1.4
    M_crit = 5 # Squared critical Mach number (above which density is modified)
    c_ref = 1
    dynP = 0.5*100 # dynamic pressure

    # mesh an airfoil
    pars = {'xLength' : 5, 'yLength' : 5, 'sEleFar' : 0.5, 'sEleAirfTE' : 0.01, 'sEleAirfLE' : 0.01}
    msh = gmsh.MeshLoader("models/n0012.geo",__file__).execute(**pars)
    pbl = f.Problem(msh, alpha, M_inf, c_ref)

    # add a medium "air"
    if M_inf == 0:
        pbl.set(f.Medium(msh, "internalField", f.Fun0EleRhoL(), f.Fun0EleMachL(), f.Fun0EleCpL(), f.Fun0PosPhiInf(alpha)))
        fCp = f.Fun0EleCpL()
    else:
        pbl.set(f.Medium(msh, "internalField", f.Fun0EleRho(gamma, M_inf, M_crit), f.Fun0EleMach(gamma, M_inf), f.Fun0EleCp(gamma, M_inf), f.Fun0PosPhiInf(alpha)))
        fCp = f.Fun0EleCp(gamma, M_inf)

    # add initial condition
    pbl.add(f.Assign(msh,"internalField", f.Fun0PosPhiInf(alpha)),"IC")
    # add boundary conditions
    pbl.add(f.Neumann(msh, "upstream", tbox.Fct1C(-U_inf[0], -U_inf[1], 0.)))
    pbl.add(f.Neumann(msh, "sideUp", tbox.Fct1C(-U_inf[0], -U_inf[1], 0.)))
    pbl.add(f.Neumann(msh, "sideLw", tbox.Fct1C(-U_inf[0], -U_inf[1], 0.)))
    pbl.add(f.Neumann(msh, "downstream", tbox.Fct1C(-U_inf[0], -U_inf[1], 0.)))
    pbl.add(f.Neumann(msh, "airfoil", tbox.Fct1C(0., 0., 0.)))
    # add Wake boundary and Kutta conditions
    pbl.add(f.Wake(msh, ["wakeUp", "wakeLw", "internalField"]))
    pbl.add(f.Kutta(msh, ["teUp", "teLw", "airfoil", "internalField"]))
    # identify the f/s boundary
    airfoil = f.Boundary(msh, "airfoil")
    pbl.add(airfoil)

    # initialize mesh deformation handler
    mshDef = tbox.MshDeformation(msh, "internalField", ["upstream", "sideUp", "sideLw", "downstream"], ["airfoil"], ["wakeUp", "wakeLw"])

    # initialize solver
    solver = f.Solver(pbl)
    initSolver(solver)
    solver.initialize()

    return Module(msh, mshDef, airfoil, solver, fCp, dynP)

    
