#!/usr/bin/env python
# -*- coding: latin-1; -*-
# Battery FSI
# Revisited for FSI by Marco Lucio CERQUAGLIA and David THOMAS
# Originally written for "Metalub" by Romain Boman
# Modified for "Metalub" by Yves CARRETTA 

import sys, glob, os, subprocess, platform

defArgs = [ r"apps".replace('/',os.sep)]
lastDir = None

def printDir(donfile):
    global lastDir
    if lastDir!=os.path.dirname(donfile):
        lastDir=os.path.dirname(donfile)

def runOne(donfile, usegui):
    global exeFile   
    logfile = donfile.replace('.py','.log')
    resfile=donfile.replace('.py','.res')
    
    #log file
    if not os.path.isfile(logfile) \
        or os.path.getmtime(logfile) < os.path.getmtime(donfile):
        printDir(donfile)
        exe='FSI'
        precise=''
        cmd = [r"python"]
        cmd += [donfile]
        cmd += [' True ']
        cmd += [usegui]
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
    if not os.path.isfile(resfile) \
        or os.path.getmtime(resfile) < os.path.getmtime(donfile):
        fres = open(resfile,'w')
        
        try:
            if os.path.isfile(logfile):
                for line in open(logfile,'r'):
                    for exp in [ 'RES-FSI-', '[Successful Run FSI]', '[cpu FSI]', '[Time steps FSI]', '[Mean n. of FSI Iterations]']:
                        if line.find(exp)!=-1:
                            fres.write("%s" % (line))        
        except: # ctrl-c
            fres.close()
            cleanOne(donfile)
            raise

        fres.close()

        if not verifOne(donfile): # check for results
            print '\tFAILURE!'
            os.utime(donfile, None) # touch donfile

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
            for exp in [ 'RES-FSI-', '[Successful Run FSI]', '[cpu FSI]', '[Time steps FSI]', '[Mean n. of FSI Iterations]']:
                if line.find(exp)!=-1:
                    tsc.append(line)
    return tsc

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

def process(args, cmd, usegui):
    global defArgs
    if not args: args=defArgs
    for donfile in loopOn(args):
        if cmd=='run':
            runOne(donfile, usegui)
        elif cmd=='clean':
            cleanOne(donfile)
        elif cmd=='rerun':
            cleanOne(donfile)
            runOne(donfile, usegui)


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

def getFSIPath():
    batteryPath = os.path.abspath(os.path.dirname(sys.argv[0]))
    import fsi_pyutils
    fsi_pyutils.findbins()
    import FSICoupler 
    pth = os.path.dirname(fsi.__file__)  
    return pth

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

usage: %s [run|rerun|verif|clean] files

examples:
  %s run: start/continue the battery
  %s rerun: restart the battery (clean before run)
  %s verif: create the summary file gathering the results
  %s clean: clean all results

""" % (exe, exe, exe, exe, exe, exe)
    print txt


def main():
    global fsiExe
    #execfile('.pythonrc.py')
    i=1
    while i<len(sys.argv):
        # "Initialisation" 
        # On définit l'extension des fichiers qui seront lus et l'executable à utiliser
        cmd=sys.argv[i]
        usegui = sys.argv[i+1]
        fsiExe = r"Python"         
        # Execution de la batterie                
        if cmd in ['run', 'clean', 'rerun']:
                process(sys.argv[i+1:], cmd, usegui)
                break
        elif cmd=='verif':
            verif(sys.argv[i+1:])
            break
        else:
            usage(os.path.basename(sys.argv[0]))
            break
        i+=1
    else:
        usage(os.path.basename(sys.argv[0]))

if __name__=="__main__":
    main()

    
