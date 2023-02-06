#! /usr/bin/env python
# -*- coding: utf-8 -*-

''' 

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

testing.py
Test classes to use CTest in CUPyDO.
Authors : David THOMAS, Marco Lucio CERQUAGLIA, Romain BOMAN, Adrien CROVATO

'''

class CTest(object):
    def __init__(self, name, val, expected, maxdiff=1e-10, forceabs=False):
        self.name     = name                  # name of the test
        self.val      = float(val)            # calculated value
        self.expected = float(expected)       # expected value
        self.maxdiff  = abs(float(maxdiff))   # tolerance on max difference

        self.forceabs = forceabs              # force the calculation of an absolute diff
        
    def run(self):
        ok = True
        adiff = abs(self.val-self.expected) # absolute diff
        
        if not self.forceabs:
            diff = adiff/abs(self.expected) # relative diff
            typ='rel'
            percent = '%3.1f%%' % (self.maxdiff*100)
        else:
            diff = adiff # absolute diff
            typ='abs'
            percent = '%f' % self.maxdiff            
 
        print("[CTest] %s = %f (expected %f +/- %s)" % \
            (self.name, self.val, self.expected, percent)) 
                   
        if diff<=self.maxdiff:
            sgn = '<='
            info = "ok"
        else:
            sgn = '>'
            ok = False
            info = "wrong!"
        print("\t%s diff = %e %s %e [%s]" % (typ, diff, sgn, self.maxdiff, info))
        return ok   
            
class CTests(object):
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
            print('All tests are OK!')
        else:
            print('\n\n')
            raise Exception('Some tests failed!')
          
    
    
