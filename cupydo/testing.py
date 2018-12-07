# -*- coding: utf-8 -*-
# testing classes (for ctest)
# Authors Romain Boman

class ccolors:
    ANSI_RED    = '\x1b[1;31m'
    ANSI_GREEN  = '\x1b[1;32m'
    ANSI_BLUE   = '\x1b[1;34m'
    ANSI_RESET  = '\x1b[0m'

class CTest:
    def __init__(self, name, val, expected, maxdiff=1e-10, refval=0.0, forceabs=False):
        self.name     = name                  # name of the test
        self.val      = float(val)            # calculated value
        self.expected = float(expected)       # expected value
        self.maxdiff  = abs(float(maxdiff))   # tolerance on max difference
        self.refval   = abs(float(refval))    # optional value used as denominator 
                                              #   if the expected one is close to 0
        self.forceabs = forceabs              # force the calculation of an absolute diff
        
    def run(self):
        ok = True
        adiff = abs(self.val-self.expected) # absolute diff
        
        # the ref value is the largest among the expected value and
        # the one provided by the user
        denom = max(abs(self.refval), abs(self.expected))
        
        if not self.forceabs and denom>self.maxdiff:
            diff = adiff/denom # relative diff
            typ='rel'
            percent = '%3.1f%%' % (self.maxdiff*100)
        else:
            diff = adiff # absolute diff
            typ='abs'
            percent = '%f' % self.maxdiff            
 
        print "[CTest] %s = %f (expected %f +/- %s)" % \
            (self.name, self.val, self.expected, percent) 
                   
        if diff<=self.maxdiff:
            sgn = '<='
            info = "ok"
        else:
            sgn = '>'
            ok = False
            info = "wrong!"
        print "\t%s diff = %e %s %e [%s]" % (typ, diff, sgn, self.maxdiff, info)
        return ok   
            
class CTests:
    def __init__(self):
        self.tests = []
    def add(self, t):
        self.tests.append(t)
    def run(self):
        ok = True
        for i,t in enumerate(self.tests):
            #print "running test %d/%d..." % (i+1,len(self.tests))
            ok = t.run() and ok
            
        if ok:
            print(ccolors.ANSI_GREEN + '** All tests are OK!' + ccolors.ANSI_RESET)
        else:
            print('\n\n')
            raise Exception(ccolors.ANSI_RED + '** Some tests failed!' + ccolors.ANSI_RESET)
          
    
    
