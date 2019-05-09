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
    def __init__(self, _msh, _mshDef, _mshWriter, _boundary, _solver, _fCp, _dynP):
        # objects
        self.msh = _msh
        self.mshDef = _mshDef
        self.mshWriter = _mshWriter
        self.boundary = _boundary
        self.solver = _solver
        self.fCp = _fCp
        # parameters
        self.dynP = _dynP

def initSolver(_solver):
    _solver.relTol = 1e-6
    _solver.absTol = 1e-8
    _solver.lsTol = 1e-6
    _solver.maxIt = 25
    _solver.maxLsIt = 10
    _solver.avThrsh = 1e-2

def getFlow():
    # define flow variables
    alpha = 3*math.pi/180
    U_inf = [math.cos(alpha), math.sin(alpha)] # norm should be 1
    M_inf = 0.3
    gamma = 1.4
    M_crit = 5 # Squared critical Mach number (above which density is modified)
    c_ref = 1
    dynP = 0.5*102*102 # dynamic pressure
    dim = len(U_inf)

    # mesh an airfoil
    pars = {'xLgt' : 5, 'yLgt' : 5, 'msF' : 1., 'msA' : 0.01}
    msh = gmsh.MeshLoader("models/diamond_fluid.geo",__file__).execute(**pars)
    gmshWriter = tbox.GmshExport(msh)
    mshCrck = tbox.MshCrack(msh, dim, gmshWriter, "wake", ["field", "field_", "airfoil", "airfoil_", "downstream", "downstream_"])
    pbl = f.Problem(msh, dim, alpha, M_inf, c_ref, c_ref, 0., 0., 0.)

    # add a medium "air"
    if M_inf == 0:
        fCp = f.Fun0EleCpL()
        pbl.set(f.Medium(msh, "field", f.Fun0EleRhoL(), f.Fun0EleMachL(), fCp, f.Fun0PosPhiInf(dim, alpha)))
    else:
        fCp = f.Fun0EleCp(gamma, M_inf)
        pbl.set(f.Medium(msh, "field", f.Fun0EleRho(gamma, M_inf, M_crit), f.Fun0EleMach(gamma, M_inf), fCp, f.Fun0PosPhiInf(dim, alpha)))

    # add initial condition
    pbl.add(f.Assign(msh,"field", f.Fun0PosPhiInf(dim, alpha)),"IC")
    # add boundary conditions
    pbl.add(f.Neumann(msh, "upstream", tbox.Fct1C(-U_inf[0], -U_inf[1], 0.)))
    pbl.add(f.Neumann(msh, "side", tbox.Fct1C(-U_inf[0], -U_inf[1], 0.)))
    pbl.add(f.Neumann(msh, "downstream", tbox.Fct1C(-U_inf[0], -U_inf[1], 0.)))
    pbl.add(f.Neumann(msh, "airfoil", tbox.Fct1C(0., 0., 0.)))
    # add wake/Kutta boundary conditions
    pbl.add(f.Wake(msh, ["wake", "wake_", "field"]))
    pbl.add(f.Kutta(msh, ["te", "wake_", "airfoil", "field"]))
    # identify the boundary
    airfoil = f.Boundary(msh, ["airfoil", "field"])
    pbl.add(airfoil)

    # initialize mesh deformation handler
    mshDef = tbox.MshDeform(msh, dim)
    mshDef.setField("field")
    mshDef.setFixed(["upstream", "side", "downstream"])
    mshDef.setMoving(["airfoil"])
    mshDef.setInternal(["wake", "wake_"])

    # initialize solver
    solver = f.Newton(pbl)
    initSolver(solver)

    return Module(msh, mshDef, gmshWriter, airfoil, solver, fCp, dynP)
