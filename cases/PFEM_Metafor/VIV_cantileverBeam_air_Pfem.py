#! /usr/bin/env python
# -*- coding: latin-1; -*-
# $Id: $

import sys, os, os.path

runPath = os.path.dirname(sys.modules[__name__].__file__)
filePath = os.path.abspath(os.path.dirname(__file__))
fileName = os.path.splitext(os.path.basename(__file__))[0]

import pfem

import pfemtools as wt
import viewer as v
    
w = None

class Module:
    def __init__(self, w, msh, pbl, scheme, extManager, gui):
       self.w = w
       self.msh = msh
       self.pbl = pbl       
       self.scheme = scheme
       self.extManager = extManager
       self.gui = gui

def getPfem():
    global contactTag
    global w
    if w: return w
    w = pfem
    
    mshFile = runPath+os.sep+'VIV_cantileverBeam.msh'
    
    rho0 = 1.18
    mu = 1.82e-5
    U = 0.513
    
    pbl = w.Problem()
    pbl.rho0 = rho0
    pbl.mu = mu
    pbl.nonLinAlgorithm = 1
    pbl.solScheme = 1
    pbl.alpha = 100.0
    pbl.extP = 0.
    pbl.scalingU = U
    
    msh = w.MshData(pbl)
    msh.load(mshFile)
    print msh
    
    contactTag = w.Tag(100, "Contact" , 2)
    msh.ptags[100] = contactTag
    msh.ntags["Contact"] = contactTag
    
    scheme = w.BackwardEuler(msh, pbl)
    
    w.Medium(msh, "Contact", 0., 0., 0., 4)
    w.Medium(msh, 23, mu, rho0, 3)
    w.Medium(msh, 22, mu, rho0, 3)
    w.Medium(msh, 25, mu, rho0, 1)
    
    print msh.media.size()
    
    # boundaries
    # w.Boundary(msh, 19, 1, 0.0)
    # w.Boundary(msh, 19, 2, 0.0)
    w.Boundary(msh, 19, 7, 0.0)
    w.Boundary(msh, 20, 3, 0.0)
    w.Boundary(msh, 20, 6, 0.0)
    w.Boundary(msh, 21, 1, U)
    w.Boundary(msh, 21, 2, 0.)
    w.Boundary(msh, 21, 5, 0.0)
    w.Boundary(msh, 22, 1, 0.0)
    w.Boundary(msh, 22, 2, 0.0)
    w.Boundary(msh, 23, 1, 0.0)
    w.Boundary(msh, 23, 2, 0.0)
    
    # --- Necessary for the correct elimination of nodes on the Outlet using slip conditions ---
    n2 = msh.getNode(2)
    n3 = msh.getNode(3)
    n2.isOnAxisOfSymmetry = True
    n3.isOnAxisOfSymmetry = True
    # ---
    
    scheme.savefreq=100
    scheme.nthreads=3
    scheme.gamma = 0.7
    scheme.omega = 0.7
    scheme.addRemoveNodesOption = True
    
    #Results
    extManager = w.ExtractorsManager(msh)
    
    gui = v.MeshViewer(msh, scheme, True) 
    
    return Module(w, msh, pbl, scheme, extManager, gui)