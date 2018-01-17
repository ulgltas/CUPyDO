#!/usr/bin/env python
# -*- coding: latin-1; -*-

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

Battery CUPyDO
Revisited for CUPyDO by Marco Lucio CERQUAGLIA and David THOMAS
Originally written for "Metalub" by Romain Boman
Modified for "Metalub" by Yves CARRETTA

'''

import sys, glob, os, subprocess, platform, argparse, time

defArgs = [ r"tests".replace('/',os.sep)]
lastDir = None

def printDir(donfile):
    global lastDir
    if lastDir!=os.path.dirname(donfile):
        lastDir=os.path.dirname(donfile)

def runOne(donfile, nbProcs):
    global exeFile   
    logfile = donfile.replace('.py','.log')
    resfile=donfile.replace('.py','.res')
    
    #log file
    if not os.path.isfile(logfile) \
        or os.path.getmtime(logfile) < os.path.getmtime(donfile):
        printDir(donfile)
        exe='FSI'
        precise=''
        if int(nbProcs) == 0:
            cmd = [r"python"]
            cmd += [donfile]
        else:
            cmd = [r"mpirun"]
            cmd += ["--np"]
            cmd += [nbProcs]
            cmd += ["python"]
            cmd += [donfile]
        
        print '\t[%s] %s => %s %s' % (exe, os.path.basename(donfile), os.path.basename(logfile), precise)
        flog = open(logfile,'w')
        
        try:
            p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=flog, stderr=flog, env=os.environ, shell=False)
            pin = p.stdin
            pin.write("__file__=r'%s'\n" % donfile )
            pin.write("execfile(r'%s')\n" % donfile)
            pin.close()
            retcode = p.wait()        
        except: # ctrl-c
            flog.close()
            cleanOne(donfile)
            raise

        flog.close()
    
        #res file
        fres = open(resfile,'w')
        
        try:
            if os.path.isfile(logfile):
                for line in open(logfile,'r'):
                    for exp in [ 'RES-FSI-', '[Successful Run FSI]', '[cpu FSI total]', '[Time steps FSI]', '[Mean n. of FSI Iterations]']:
                        if line.find(exp)!=-1:
                            fres.write("%s" % (line))        
        except: # ctrl-c
            fres.close()
            cleanOne(donfile)
            raise
        
        fres.close()
        
        if not verifOne(donfile): # tsc or not checkOneRun(tsc): # check for results
            print '\tFAILURE!'
            os.utime(donfile, (time.time()+1.0,time.time()+1.0)) # touch donfile

def cleanOne(donfile):
    for ext in ['.log','py.log','.res','.err','.pyc']:
        rmfile = donfile.replace('.py',ext)
        if os.path.isfile(rmfile):
            printDir(donfile)
            print "\trm %s" % os.path.basename(rmfile)
            os.remove(rmfile)

def verifOne(donfile):
    tsc=[]
    resfile=donfile.replace('.py','.res')
    if os.path.isfile(resfile):
        for line in open(resfile,'r'):
            for exp in [ 'RES-FSI-', '[Successful Run FSI]', '[cpu FSI total]', '[Time steps FSI]', '[Mean n. of FSI Iterations]']:
                if line.find(exp)!=-1:
                    tsc.append(line)
    
    if not checkOneRun(tsc):
        tsc=[]
    
    return tsc

def checkOneRun(tsc):
    runOk = False
    if any('[Successful Run FSI]: True' in s for s in tsc):
        runOk = True
    return runOk

def loopOnOne(file):
    if os.path.isdir(file):
        if os.path.basename(file) == '.svn':
            return
        subfiles = os.listdir(file)
        subfiles.sort()
        for sfile in subfiles:
            for f in loopOnOne(os.path.join(file, sfile)):
                yield f
    elif os.path.isfile(file):
        if os.path.basename(file) == '__init__.py':
            return
        elif os.path.splitext(file)[1]=='.py':
            yield file

def loopOn(files):
    if not files: return
    files.sort()
    for wfile in files:
        for file in glob.glob(wfile):
            for f in loopOnOne(file):
                yield f

def process(args, cmd, nbProcs):
    global defArgs
    if not args: args=defArgs
    for donfile in loopOn(args):
        if cmd=='run':
            runOne(donfile, nbProcs)
        elif cmd=='clean':
            cleanOne(donfile)
        elif cmd=='rerun':
            cleanOne(donfile)
            runOne(donfile, nbProcs)


def machineid():
    uname = platform.uname()
    if uname[0] == 'Windows' or uname[2] == 'Windows' or uname[0].find("CYGWIN")!=-1:
        return "Windows"
    elif uname[0] == 'Linux':
        if uname[4] == 'x86_64':
            return "Linux64"
        elif uname[4] == 'alpha':
            return "AlphaLinux"
        else:
            return "Linux"
    else:
        return uname[0]

def verif(args):
    print "gathering the results in fsi/verif"
    global defArgs
    if not args: args=defArgs
    
    fsiPth = ('..')
    verifPth = ('..') + os.sep + 'verif'
    
    files={}
    fext = ['failed', 'results']
    for ext in fext:
        fname = ( verifPth + os.sep + "%s-%s.ascii" % (ext, machineid())).replace('/',os.sep)               
        files[ext] = open(fname,'w')

    ntest=0
    nfailed=0
    for donfile in loopOn(args):
        ntest+=1
        tsc = verifOne(donfile)
        if not tsc:
            files['failed'].write("%s: FAILED!\n" % donfile.replace(fsiPth,''))
            nfailed+=1
        for line in tsc:
            files['results'].write("%s: %s" %(donfile.replace(fsiPth,''), line))

    print "%d/%d tests OK (%d failed)" %(ntest-nfailed, ntest, nfailed)
    for ext in fext:
        files[ext].close()

def usage(exe):
    txt="""
%s : Battery script for FSI

usage: %s [run|rerun|verif|clean] tests

examples:
  %s run: start/continue the battery
  %s rerun: restart the battery (clean before run)
  %s verif: create the summary file gathering the results
  %s clean: clean all results

""" % (exe, exe, exe, exe, exe, exe)
    return txt


def main():
    
    # Analyzing command line options
    parser = argparse.ArgumentParser(usage=usage(os.path.basename(sys.argv[0])))
    parser.add_argument('cmd', type=str, choices=['run','rerun','clean','verif'],
                        help='instruction to be executed by battery')
    parser.add_argument('--np', dest='nbProcs', type=str, default='0',
                        help='number of processors for MPI computations')
    parser.add_argument('tests', type=str, nargs='+',
                        help='test files name')
    args = parser.parse_args()
    
    cmd = args.cmd
    nbProcs = args.nbProcs
    tests = args.tests
    
    # Execution of the battery                
    if cmd in ['run', 'clean', 'rerun']:
        process(tests, cmd, nbProcs)
    elif cmd=='verif':
        verif(tests)
    else:
        parser.print_help()

if __name__=="__main__":
    main()

    
