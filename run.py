#!/usr/bin/env python
# -*- coding: utf8 -*-
# test encoding: à-é-è-ô-ï-€

# Run script for CUPyDO
# Romain Boman and Adrien Crovato

import os
import sys
import json
import cupydo.utilities as cupyutil

# These two classes will not redirect c++ output without the underlying c++ code (waves/fwk)
class DupStream(object):
    def __init__(self, stream1, stream2):
        self.stream1 = stream1
        self.stream2 = stream2

    def write(self, data):
        self.stream1.write(data)
        self.stream2.write(data)

    def flush(self):
        self.stream1.flush()
        self.stream2.flush()

class Tee(object):
    def __init__(self, name):
        self.file = open(name, 'w')
        self.stdoutbak = sys.stdout
        self.stderrbak = sys.stderr
        sys.stdout = DupStream(sys.stdout, self.file)
        sys.stderr = DupStream(sys.stderr, self.file)

    def __del__(self):
        sys.stdout = self.stdoutbak
        sys.stderr = self.stderrbak
        self.file.close()

def addPath(p):
    # windows binaries are located somewhere else
    pathw = os.path.join(p, 'Release')
    if os.path.isdir(pathw):
        p = pathw
    # add folder to path if it exists
    if os.path.isdir(p):
        print('INFO: adding %s to PYTHONPATH' % p)
        sys.path.append(p)
    else:
        print('INFO: %s not found!' % p)

def setPath():
    # Set paths
    cupdir = os.path.abspath(os.path.split(__file__)[0])
    topdir = os.path.abspath(os.path.dirname(cupdir))
    haveMPI, comm, myid, numberPart = cupyutil.getMpi()

    try: file = open('solverPaths.json')
    except: raise Exception('Could not find solverPaths.json')
    solverPaths = json.load(file)

    try: addPath(solverPaths['Metafor'])
    except: pass
    try: addPath(solverPaths['PFEM3D'])
    except: pass
    try: addPath(solverPaths['Modali'])
    except: pass
    try: addPath(solverPaths['VLM'])
    except: pass
    try: addPath(solverPaths['Waves'])
    except: pass
    try: addPath(solverPaths['PFEM'])
    except: pass
    try: addPath(solverPaths['SU2'])
    except: pass
    try: addPath(solverPaths['NativeSolid'])
    except: pass
    try: addPath(solverPaths['pyBeam'])
    except: pass

    if myid == 0:
        print('PYTHONPATH = %s\n' % sys.path)
    cupyutil.mpiBarrier(comm)

def main():
    # Global variables
    global __file__

    # Find solvers and set paths
    setPath()

    # Parse arguments
    args = cupyutil.parseArgs()

    # MPI comm to avoid repetition
    haveMPI, comm, myid, numberPart = cupyutil.getMpi()
    # Process
    for file in args.file:
        if not os.path.isfile(file):
            raise RuntimeError('Could not find file:', file)
        else:
            file = os.path.abspath(file)
            # change dir
            cupyutil.setDirs(file)
            # split streams
            __file__ = file
            cupyutil.mpiPrint("[run.py] __file__" + __file__, comm)
            tee = Tee('stdout.txt')
            # run
            import time, platform
            cupyutil.mpiPrint('-' * 80, comm)
            cupyutil.mpiPrint('Starting test:' + file, comm)
            cupyutil.mpiPrint('Time:' + time.strftime('%c'), comm)
            cupyutil.mpiPrint('Hostname:' + platform.node(), comm)
            exec(open(file, 'r', encoding='utf8').read(), globals())

if __name__ == '__main__':
    main()