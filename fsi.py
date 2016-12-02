#!/usr/bin/env python
# -*- coding: latin-1; -*-
# ~/dev/Metafor/oo_metaB/bin/Metafor -nogui ./fsi.py

import math
import numpy as np
import os, os.path, sys, time, string

import socket, fnmatch
import pyutils

# global vars (underscore prevent them to be imported with "from module import *")
_theModule  = None
_theWDir    = None # workspace directory
_theWDirRoot = os.getcwd()  # base directory du workspace

# ------------------------------------------------------------------------------

def load(fsiTxt):
    """ load a module and make it the current one
    """
    global _theModule
    
    if not isinstance(fsiTxt, str): raise Exception("argument must be a string!")
    if _theModule: raise Exception("no more than one module may be loaded!")

    if fsiTxt=="__main__": # on est dans un script => getpfem est/sera dans __main__
        _theModule = sys.modules['__main__']
    else: 
        if os.path.isfile(fsiTxt): # c'est une nom de fichier => on convertit en nom de module
            module = pyutils.fileToModule(fsiTxt, verb=False)
            if not module: raise Exception('file "%s" is not a reachable python module!' % fsiTxt)
            fsiTxt = module
        # ici, on a un nom de module...
        exec("import %s" % fsiTxt) # peut on le faire plus tard?
        exec("globals()['_theModule'] = %s" % fsiTxt)
        print "module '%s' loaded!" % fsiTxt
        #print '_theModule', _theModule
        
    setTheWDir(fsiTxt)
    setDir(_theWDir)
    _chkWorkspace()
    
    return _theModule
    
def setTheWDir(fsiTxt): 
    global _theWDir
    _theWDir = defaultWDir(fsiTxt)
    # on fait un chdir si ce rep existe!
    if os.path.isdir(_theWDir):
        setDir(_theWDir)
        
    return _theWDir
     
def defaultWDir(moduleTxt):
    global _theWDirRoot
    return os.path.join(_theWDirRoot, os.path.join("workspace", moduleTxt.replace('.','_') )) # default WDir for the module
    
def setDir(wdir):  # (avant instance) - change le rep courant
    global _theWDirRoot
    """ change the current working directory
    """
    if wdir:
        if not os.path.isabs(wdir):
            os.path.join(_theWDirRoot, wdir)
        try:
            if not os.path.isdir(wdir):
                os.makedirs(wdir)
            os.chdir(wdir)
        except OSError, e:
            print "OSError : ", e
            # check UAC (Vista) or root-protected dirs (Unix)    
            text = 'Directory "%s" may not be created.\n' % wdir
            text += 'Disable UAC (vista) or ask for root privilege (Unix) if you still want to work here.\n'
            text += 'Otherwise, rebase to another directory!' 
            try: wrap.GUILink().warningMsg(text)
            except: pass
            raise Exception('directory "%s" cannot be created' % wdir)    
    global _theWDir
    _theWDir = os.getcwd() # stoque un chemin absolu meme si ca a merdé.
    print "_theWDir = ", _theWDir
    
def _chkWorkspace():
    """ Check/Clean workspace
    """
    flist=[]
    for file in os.listdir('.'):
        if os.path.isdir(file):
            raise Exception("Your workspace contains one or more directories!")
        elif os.path.isfile(file):
            for sk in ['*.fdb','*.conf','*.msh']: # files to keep
                if fnmatch.fnmatch(os.path.basename(file), sk):
                    break
            else:
                flist.append(os.path.abspath(file))
    
    if flist:
        answer = True
        try: answer = wrap.GUILink().askDelete(flist)
        except: pass
        if answer:
            for file in flist:
                try:
                    os.remove(file)
                except:
                    print "Unable to remove", file
        else:
            raise Exception("meta() cancelled!")

# ----------------------------------------------------------------------
#  Generic solid solver class
# ----------------------------------------------------------------------

class SolidSolver:
    """
    Des.
    """

    def __init__(self):

        self.haloNodeList = {}

        # --- Create the array for external communication (displacement, velocity and velocity at the previous time step --- #
        self.nodalDisp_X = np.zeros(self.nPhysicalNodes)
        self.nodalDisp_Y = np.zeros(self.nPhysicalNodes)
        self.nodalDisp_Z = np.zeros(self.nPhysicalNodes)
        self.nodalVel_X = np.zeros(self.nPhysicalNodes)
        self.nodalVel_Y = np.zeros(self.nPhysicalNodes)
        self.nodalVel_Z = np.zeros(self.nPhysicalNodes)
        self.nodalVel_XNm1 = np.zeros(self.nPhysicalNodes)
        self.nodalVel_YNm1 = np.zeros(self.nPhysicalNodes)
        self.nodalVel_ZNm1 = np.zeros(self.nPhysicalNodes)

    def setInitialDisplacements(self):
        return

    def preprocessTimeIter(self, timeIter):
        return

    def run(self):
        return

    def __setCurrentState(self):
        return

    def getNodalDisplacements(self):
        """
        Des.
        """
        
        return (self.nodalDisp_X, self.nodalDisp_Y, self.nodalDisp_Z)

    def getNodalInitialPositions(self):
        return


    def getNodalVelocity(self):
        """
        des.
        """
        
        return (self.nodalVel_X, self.nodalVel_Y, self.nodalVel_Z)
    
    def getNodalVelocityNm1(self):
        """
        Des.
        """
        
        return (self.nodalVel_XNm1, self.nodalVel_YNm1, self.nodalVel_ZNm1)
    
    def getNodalIndex(self, iVertex):
        return
    
    def fakeFluidSolver(self, time):
        return
    
    def applyNodalLoads(self, load_X, load_Y, load_Z, time):
        return
    
    def update(self):

        self.nodalVel_XNm1 = self.nodalVel_X.copy()
        self.nodalVel_YNm1 = self.nodalVel_Y.copy()
        self.nodalVel_ZNm1 = self.nodalVel_Z.copy()
    
    def save(self):
        return
    
    def initRealTimeData(self):
        return
    
    def saveRealTimeData(self, time, nFSIIter):
        return

    def remeshing(self):
        return

    def exit(self):
        return

# ----------------------------------------------------------------------
#  Generic fluid solver class
# ----------------------------------------------------------------------

class FluidSolver:
    """
    Des.
    """

    def __init__(self):

        self.haloNodeList = {}

        self.nodalLoad_X = np.zeros((self.nPhysicalNodes))
        self.nodalLoad_Y = np.zeros((self.nPhysicalNodes))
        self.nodalLoad_Z = np.zeros((self.nPhysicalNodes))

    def setInitialMeshDeformation(self):
        return

    def preprocessTimeIter(self, timeIter):
        return
    
    def run(self):
        return
    
    def getNodalIndex(self, iVertex):
        return
    
    def fakeSolidSolver(self, time):
        return
    
    def getNodalLoads(self):

        return (self.nodalLoad_X, self.nodalLoad_Y, self.nodalLoad_Z)
    
    def applyNodalDisplacements(self, dx, dy, dz, dx_nM1, dy_nM1, dz_nM1, haloNodesDisplacements,time):
        return

    def update(self, dt):
        return
    
    def save(self, nt):
        return
    
    def initRealTimeData(self):
        return
    
    def saveRealTimeData(self, time, nFSIIter):
        return
    
    def remeshing(self):
        return

    def meshUpdate(self, nt):
        return

    def exit(self):
        return

# ----------------------------------------------------------------------
# Convergence criteria for implicit FSI coupling 
# ----------------------------------------------------------------------

class FsiCriterion:
    def __init__(self, toll):
        self.error = 1.0
        self.toll = toll
    
    def isVerified(self):
        if self.error < self.toll:
            return True
        else:
            return False
            
class SolidPosErrorCriterion(FsiCriterion):
    def __init__(self, solidNodes, toll):
        FsiCriterion.__init__(self, toll)
        self.pos = {}
        for no in solidNodes.iterkeys():
            node = solidNodes[no]
            px = node.pos.x[0]
            py = node.pos.x[1]
            self.pos[no] = (px, py)
            '''node,fx,fy = solidNodes[no]
            px = node.getPos(Configuration().currentConf).get1() # current x          
            py = node.getPos(Configuration().currentConf).get2() # current y
            self.pos[no] = (px, py)'''
            
    def update(self, solidNodes):
        err2 = 0.
        pos2 = 0.
        for no in solidNodes.iterkeys():
            node,fx,fy = solidNodes[no]
            px = node.getPos(Configuration().currentConf).get1() # current x          
            py = node.getPos(Configuration().currentConf).get2() # current y
            pxOld = self.pos[no][0]
            pyOld = self.pos[no][1]
            
            err2 += (px - pxOld)*(px - pxOld) + (py - pyOld)*(py - pyOld)
            pos2 += pxOld*pxOld + pyOld*pyOld
            
            self.pos[no] = (px, py)
        self.error = math.sqrt(err2)/math.sqrt(pos2)

class DispResidualBasedCriterion(FsiCriterion):
    def __init__(self, toll):
        FsiCriterion.__init__(self, toll)
        
    def update(self, res):
        neq = len(res)
        self.error = (1.0/(math.sqrt(neq)))*(np.linalg.norm(res))

# ----------------------------------------------------------------------
# Partitioned FSI coupling algorithms
# ----------------------------------------------------------------------

# --- Generic FSI algorithm ---

class FsiAlgorithm:
    def __init__(self, _solid, _fluid, _dt, _tTot, _criterion):
        self.solid = _solid
        self.fluid = _fluid
        self.dt = _dt
        self.tTot = _tTot
        self.fsiCriterion = _criterion
        
        self.t = 0.
    
    def run(self, t1, t2):
        print '\nERROR: call to the generic class fsiAlgorithm run function!\n'
        exit(1)
        
    def remeshing(self):
        self.solid.remeshing()
        self.fluid.remeshing()


# --- Fixed-point algorithm with Aitken's relaxation ---

class FixedPointAitkenRelaxationAlgorithm(FsiAlgorithm):
    def __init__(self, _solid, _fluid, _dt, _tTot, _criterion, omega_m):
        FsiAlgorithm.__init__(self, _solid, _fluid, _dt, _tTot, _criterion)
        self.omega_max = omega_m
        
    def run(self):
    
        # Definition of necessary quantities for relaxation
        omega = self.omega_max
        omega0 = self.omega_max
        omega_n = self.omega_max
        residual = np.zeros(3*len(self.solid.fnods))
        residual0 = np.zeros(3*len(self.solid.fnods))
        d = np.zeros(3*len(self.solid.fnods))
        d0 = np.zeros(3*len(self.solid.fnods))
        dn = np.zeros(3*len(self.solid.fnods))
        dn_1 = np.zeros(3*len(self.solid.fnods))
        dn_2 = np.zeros(3*len(self.solid.fnods))
        d_relaxed = np.zeros(3*len(self.solid.fnods))
        d_guess = np.zeros(3*len(self.solid.fnods))
        
        t1 = 0.0  # initial time
        t2 = t1
        nt = 0 # number of time steps
        
        logFile = open('fsi.log', 'w')        
        
        while self.t < self.tTot: 
            logFile.write('\n')
            logFile.write('--- New time step: nt=%d, time=%f ---\n' % (nt, self.t))
            
            logFile.write('Initializing FSI coupling...\n')
            
            # Solid predictor
            if nt == 0: 
                d_guess = dn # Order 0 prediction
            elif nt == 1:
                d_guess = 2.0*dn - dn_1 # Order 1 prediction
            else:
                d_guess = (5.0/2.0)*dn - 2.0*dn_1 + 0.5*dn_2 # Order 2 prediction
                
            # Fluid solver step
            dx, dy, dz = self.solid.getNodalDisplacements()
            dNm1 = np.split(dn, 3)

            self.fluid.applyNodalDisplacements(dx, dy, dz, dNm1[0], dNm1[1], dNm1[2], [] , t2)
            
            t2=t1+self.dt  # time to be calculated
            nFsiIt = 0 
            
            self.fluid.run(t1,t2)
            if not self.fluid.runOK:
                print '\nERROR: fluid solver did not converge!\n'
                exit(1)
                    
            # Solid solver step
            fx, fy, fz = self.fluid.getNodalLoads()
            self.solid.applyNodalLoads(fx,fy,fz,t2)
            self.solid.run(t1,t2)
            if not self.solid.runOK:
                print '\nERROR: solid solver did not converge!\n'
                exit(1)
            dx, dy, dz = self.solid.getNodalDisplacements()
            d = np.concatenate([dx,dy,dz])
            
            # Compute residual
            residual = d - d0 
            
            # Compute FSI error
            self.fsiCriterion.update(residual)
            
            np.copyto(residual0, residual)
            
            while not self.fsiCriterion.isVerified():

                # Computation of Aitken's relaxation parameter
                if nFsiIt == 0:
                    omega = self.omega_max # omega_max is user defined
                    # omega = math.copysign(min(math.fabs(omega_n),omega_max), omega_n) # omega_max is user defined
                else:
                    omega = -omega0*((np.dot(residual0,(residual-residual0)))/(np.linalg.norm(residual - residual0)*np.linalg.norm(residual - residual0)))
                
                # Relaxation step
                d_relaxed = omega*d + (1.0-omega)*d0
                
                nFsiIt+=1
                
                logFile.write('\tFSI coupling iteration %d' % nFsiIt)
                logFile.write('\t-\trelaxation parameter = %e' % omega)
                
                # Fluid solver step
                
                d_r = np.split(d_relaxed, 3)
                
                self.fluid.applyNodalDisplacements(d_r[0], d_r[1], d_r[2], dNm1[0], dNm1[1], dNm1[2], [] , t2)
                self.fluid.run(t1,t2)
                if not self.fluid.runOK:
                    print '\nERROR: fluid solver did not converge!\n'
                    exit(1)
                
                # Solid solver step
                fx, fy, fz = self.fluid.getNodalLoads()
                self.solid.applyNodalLoads(fx,fy,fz,t2)
                self.solid.run(t1,t2)
                if not self.solid.runOK:
                    print '\nERROR: solid solver did not converge!\n'
                    exit(1)
                dx, dy, dz = self.solid.getNodalDisplacements()
                d = np.concatenate([dx,dy,dz])
                
                # Compute residual
                np.copyto(residual0, residual)
                residual = d - d0 # or d_relaxed - d?
                
                # Compute FSI error
                self.fsiCriterion.update(residual)
                logFile.write('\t-\tFSI criterion error = %e\n' % self.fsiCriterion.error)
                
                # Update variables
                omega0 = omega
                np.copyto(d0, d)

            logFile.write('Coupling convergence reached in %d iterations' % nFsiIt)
            logFile.write('\n')
            
            omega_n = omega
            np.copyto(dn_2,dn_1)
            np.copyto(dn_1,dn)
            np.copyto(dn, d)
            
            self.solid.update()
            self.fluid.update(self.dt)
            
            self.solid.save()
            self.fluid.save(nt+1)
            
            self.remeshing()
            
            t1=t2 # fsi loop has converged - time t2 is accepted'''
            nt+=1
            self.t+=self.dt
            
            self.solid.saveRealTimeData(self.t, nFsiIt)

        # end.
        logFile.close()
        
        self.solid.exit()
        self.fluid.exit()


# --- Solid test algorithm ---

class FsiSolidTestAlgorithm:
    def __init__(self, _solid):
        self.solid = _solid
        
    def run(self):
        # --------------------------
        # fake FSI solver
        # --------------------------
        
        t1 = 0.0  # initial time
        dt = 0.5  # time step size
        
        # we want 5 time steps
        for j in range(5):
        
            # each time step is arbitrarily calculated twice (for testing purpose)
            for i in range(2):
                
                t2=t1+dt  # time to be calculated
                
                self.solid.fakefluidsolver(t2)  # creates some dummy loads for time=t2
                # must be replaced by a loop over the nodes of the interface and
                # several calls to "solid.applyload(node_no, force_x, force_y)"
                
                # run solid solver
                print '='*80
                print "running from %f to %f: try #%d" % (t1,t2,i+1)
                print '='*80
                self.solid.run(t1,t2)
                
                # gets the deformed interface
                out = self.solid.getdefo()
                print out   # <= todo: give 'out' to the fluid solver and update the fluid mesh
    
            self.solid.nextstep()
            t1=t2 # fsi loop has converged - time t2 is accepted
        
        # end.
