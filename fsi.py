#!/usr/bin/env python
# -*- coding: latin-1; -*-
# ~/dev/Metafor/oo_metaB/bin/Metafor -nogui ./fsi.py

import math
import numpy as np
import os, os.path, sys, time, string

# Defining all the necessary paths for the fluid and solid solvers to run
mtfPath = 'D:/Dev/Officiel/Metafor'
pfemPath = 'D:/PFEM/pfemB/bin/Release'
pfemToolsPath = 'D:/PFEM/pfem/tools' 

sys.path.append(mtfPath)
sys.path.append(pfemPath)
sys.path.append(pfemToolsPath)

# Importing the necessary packages
from wrap import *
import toolbox.fac as fac

import socket, fnmatch
import pyutils

# global vars (underscore prevent them to be imported with "from module import *")
_thePfem = None
_theModule  = None
_theWDir    = None
_theBaseDir = os.getcwd()  # base directory
_theWDirRoot = os.getcwd()  # base directory du workspace

# ------------------------------------------------------------------------------

def load(pfemTxt):
    """ load a module and make it the current one
    """
    global _theModule
    
    if not isinstance(pfemTxt, str): raise Exception("argument must be a string!")
    if _theModule: raise Exception("no more than one module may be loaded!")

    if pfemTxt=="__main__": # on est dans un script => getpfem est/sera dans __main__
        _theModule = sys.modules['__main__']
    else: 
        if os.path.isfile(pfemTxt): # c'est une nom de fichier => on convertit en nom de module
            module = pyutils.fileToModule(pfemTxt, verb=False)
            if not module: raise Exception('file "%s" is not a reachable python module!' % pfemTxt)
            pfemTxt = module
        # ici, on a un nom de module...
        exec("import %s" % pfemTxt) # peut on le faire plus tard?
        exec("globals()['_theModule'] = %s" % pfemTxt)
        print "module '%s' loaded!" % pfemTxt
        #print '_theModule', _theModule
        
    setTheWDir(pfemTxt)
    setDir(_theWDir)
    _chkWorkspace()
    
    return _theModule
    
def setTheWDir(pfemTxt): 
    global _theWDir
    _theWDir = defaultWDir(pfemTxt)
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

# ------------------------------------------------------------------------------

class NLoad:
    """
    Nodal load
    """
    def __init__(self, val1, t1, val2, t2):
        self.val1 = val1
        self.t1 = t1
        self.val2 = val2
        self.t2 = t2
    def __call__(self, time):
        return self.val1 + (time-self.t1)/(self.t2-self.t1)*(self.val2-self.val1)
    def nextstep(self):
        self.t1 = self.t2
        self.val1 = self.val2 
   
# ------------------------------------------------------------------------------

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
        self.interfaceDisp_X = np.zeros(self.nLocalInterfacePhysicalNodes)
        self.interfaceDisp_Y = np.zeros(self.nLocalInterfacePhysicalNodes)
        self.interfaceDisp_Z = np.zeros(self.nLocalInterfacePhysicalNodes)
        self.interfaceVel_X = np.zeros(self.nLocalInterfacePhysicalNodes)
        self.interfaceVel_Y = np.zeros(self.nLocalInterfacePhysicalNodes)
        self.interfaceVel_Z = np.zeros(self.nLocalInterfacePhysicalNodes)
        self.interfaceVel_XNm1 = np.zeros(self.nLocalInterfacePhysicalNodes)
        self.interfaceVel_YNm1 = np.zeros(self.nLocalInterfacePhysicalNodes)
        self.interfaceVel_ZNm1 = np.zeros(self.nLocalInterfacePhysicalNodes)

    def setInitialDisplacements(self):

        return 0

    def preprocessTimeIter(self, timeIter):

        return 0

    def run(self):

        return 0

    def __setCurrentState(self):

        return 0

    def getInterfaceNodalDisplacement(self):
        """
        Des.
        """

        return (self.interfaceDisp_X, self.interfaceDisp_Y, self.interfaceDisp_Z)

    def getInterfaceNodalInitialPositions(self):

        return 0

    def getInterfaceNodalVelocity(self):
        """
        des.
        """

        return (self.interfaceVel_X, self.interfaceVel_Y, self.interfaceVel_Z)

    def getInterfaceNodalVelocityNm1(self):
        """
        Des.
        """

        return (self.interfaceVel_XNm1, self.interfaceVel_YNm1, self.interfaceVel_ZNm1)

    def getInterfaceNodeIndex(self, iVertex):

        return 0

    def fakeFluidSolver(self, time):

        return 0

    def applyLoad(self, load_X, load_Y, load_Z, time):

        return 0

    def update(self):

        self.interfaceVel_XNm1 = self.interfaceVel_X.copy()
        self.interfaceVel_YNm1 = self.interfaceVel_Y.copy()
        self.interfaceVel_ZNm1 = self.interfaceVel_Z.copy()

    def save(self):

        return 0

    def initRealTimeData(self):

        return 0

    def saveRealTimeData(self, time, nFSIIter):

        return 0

    def remeshing(self):

        return 0

    def exit(self):

        return 0


# --- Pfem interface --- 

class PfemSolver:
    def __init__(self, testname, bndno):
        self.testname = testname  # string (name of the module of the fluid model)
        self.bndno = bndno        # phygroup# of the f/s interface
        
        # internal vars
        self.vnods = {}           # dict of interface nodes / prescribed velocities
        self.t1      = 0.0        # last reference time        
        self.t2      = 0.0        # last calculated time
                 
        # loads the python module
        #load(self.testname)         # use toolbox.utilities
        exec("import %s" % self.testname)
        exec("module = %s" % self.testname) # link to Pfem object

        # create an instance of Pfem
        self.pfem = module.getPfem()
        
        # retrieve the f/s boundary and the related nodes
        gr = self.pfem.w.Group(self.pfem.msh, bndno)
        
        # builds a list (dict) of interface nodes
        for e in gr.tag.elems:
            for n in e.nodes:
                no = n.no
                self.vnods[no] = (n)
        
        # Pfem scheme initialization
        self.V = self.pfem.w.DoubleVector()
        self.V0 = self.pfem.w.DoubleVector()
        self.u = self.pfem.w.DoubleVector()
        self.v = self.pfem.w.DoubleVector()
        self.p = self.pfem.w.DoubleVector()
        self.velocity = self.pfem.w.DoubleVector()
        self.pfem.scheme.init(self.V,self.V0,self.u,self.v,self.p,self.velocity)
        
        self.runOK = True
        
    def run(self, t1, t2):
        """
        calculates one increment from t1 to t2.
        """
        self.runOK = self.pfem.scheme.runOneTimeStep(self.V,self.V0,self.u,self.v,self.p,self.velocity)
        
    def fakesolidsolver(self, time):
        """
        calculate some dummy positions and velocities as a function of timestep.
        it should be replaced by the solid solver in practice.
        for each node, the fsi solver may call the "fluid.applyposition" function.
        """
        out = {}
        for no in self.vnods.iterkeys():
            node = self.vnods[no]                 
            vx = -0.5 # current vx          
            vy = 0. # current vy
            px = node.pos0.x[0] + vx*self.pfem.scheme.dt # current x         
            py = node.pos0.x[1] + vy*self.pfem.scheme.dt# current y              
            out[no] = (px,py,vx,vy)
            
        self.applydefo(out)
    
    def applydefo(self, fromsolid, t2):
        """
        prescribes given nodal positions and velocities coming from solid solver to node #no
        """
        if self.pfem.scheme.t < t2:
            self.pfem.scheme.resetNodalPositions()
        
        for no in self.vnods.iterkeys():
            node = self.vnods[no]
            node.imposedU = fromsolid[no][4]
            node.imposedV = fromsolid[no][5]
            
    def getload(self, no):
        """
        returns the updated loads of the interface nodes.
        """
        node = self.vnods[no]
        valx = -node.Fint.x[0]
        valy = -node.Fint.x[1]
        return valx, valy
    
    def update(self, dt):
        self.pfem.scheme.t+=dt
        self.pfem.scheme.nt+=1
        self.pfem.scheme.updateSolutionVectors(self.V,self.u,self.v,self.p,self.velocity)
        
    def save(self, nt):
        if nt%self.pfem.scheme.savefreq==0:
            self.pfem.scheme.archive()
        self.pfem.scheme.vizu(self.u,self.v,self.p)
        
    def remeshing(self):
        self.pfem.scheme.remeshing(self.V,self.V0,self.p)
        self.pfem.scheme.updateData(self.V0,self.V)
    
    def exitFsi(self):
        self.pfem.gui.save2vtk()


# --- Data managing tools for FSI ---

def fillDispArrayFromSolidDefo(fromsolid):
    nbNodes = len(fromsolid)
    i = 0
    dlist = []
    for no in fromsolid.iterkeys():
        dlist.insert(i, fromsolid[no][2])
        dlist.insert(i+nbNodes, fromsolid[no][3])
        i+=1
    d = np.asarray(dlist)
    return d

def getDefoFromDisplacement(nods, d, d0, dt):
    v = (1.0/dt)*(d - d0)
    out = {}
    i = 0
    nbNods = len(nods)
    for no in nods.iterkeys():
        vx = v[i]
        vy = v[i+nbNods]
        out[no] = (0., 0., 0., 0., vx, vy)
        i+=1
    return out


# --- Implicit FSI coupling convergence criteria ---

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


# --- FSI partitioned coupling algorithms ---

class fsiAlgorithm:
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


# Fixed-point algorithm with Aitken's relaxation

class fixedPointAitkenRelaxationAlgorithm(fsiAlgorithm):
    def __init__(self, _solid, _fluid, _dt, _tTot, _criterion, omega_m):
        fsiAlgorithm.__init__(self, _solid, _fluid, _dt, _tTot, _criterion)
        self.omega_max = omega_m
        
    def run(self):
    
        # Definition of necessary quantities for relaxation
        omega = self.omega_max
        omega0 = self.omega_max
        omega_n = self.omega_max
        residual = np.zeros(2*len(self.solid.fnods))
        residual0 = np.zeros(2*len(self.solid.fnods))
        d = np.zeros(2*len(self.solid.fnods))
        d0 = np.zeros(2*len(self.solid.fnods))
        dn = np.zeros(2*len(self.solid.fnods))
        dn_1 = np.zeros(2*len(self.solid.fnods))
        dn_2 = np.zeros(2*len(self.solid.fnods))
        d_relaxed = np.zeros(2*len(self.solid.fnods))
        d_guess = np.zeros(2*len(self.solid.fnods))
        
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
            defo = getDefoFromDisplacement(self.solid.fnods, d_guess, d, self.dt)
            self.fluid.applydefo(defo, t2)
            
            t2=t1+self.dt  # time to be calculated
            nFsiIt = 0 
            
            self.fluid.run(t1,t2)
            if not self.fluid.runOK:
                print '\nERROR: fluid solver did not converge!\n'
                exit(1)
                    
            # Solid solver step
            for no in self.solid.fnods.iterkeys():
                valx, valy = self.fluid.getload(no)
                self.solid.applyload(no,valx,valy,t2)
            self.solid.run(t1,t2)
            if not self.solid.runOK:
                print '\nERROR: solid solver did not converge!\n'
                exit(1)
            defo = self.solid.getdefo()
            d = fillDispArrayFromSolidDefo(defo)
            
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
                defo = getDefoFromDisplacement(self.solid.fnods, d_relaxed, dn, self.dt)
                self.fluid.applydefo(defo, t2)
                self.fluid.run(t1,t2)
                if not self.fluid.runOK:
                    print '\nERROR: fluid solver did not converge!\n'
                    exit(1)
                
                # Solid solver step
                for no in self.solid.fnods.iterkeys():
                    valx, valy = self.fluid.getload(no)
                    self.solid.applyload(no,valx,valy,t2)
                self.solid.run(t1,t2)
                if not self.solid.runOK:
                    print '\nERROR: solid solver did not converge!\n'
                    exit(1)
                defo = self.solid.getdefo()
                d = fillDispArrayFromSolidDefo(defo)
                
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
            
        # end.
        logFile.close()
        
        self.solid.exitFsi()
        self.fluid.exitFsi()
        
class fsiSolidTestAlgorithm:
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
                
                # runs solid solver
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
