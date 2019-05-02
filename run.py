#!/usr/bin/env python
# -*- coding: utf8 -*-
# test encoding: à-é-è-ô-ï-€

# Run script for CUPyDO
# External solver dir should place next to CUPyDO dir
# Romain Boman and Adrien Crovato

import os
import sys

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
        import sys
        self.file = open(name, 'w')
        self.stdoutbak = sys.stdout
        self.stderrbak = sys.stderr
        sys.stdout = DupStream(sys.stdout, self.file)
        sys.stderr = DupStream(sys.stderr, self.file)

    def __del__(self):
        import sys
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
        print 'INFO: adding %s to PYTHONPATH'
        sys.path.append(p)
    else:
        print 'INFO: %s not found!'


def setPath():
    import os, sys
    # Set paths
    cupdir = os.path.abspath(os.path.split(__file__)[0])
    topdir = os.path.abspath(os.path.dirname(cupdir))

    addPath(os.path.join(topdir, 'Metafor', 'oo_metaB', 'bin'))
    addPath(os.path.join(topdir, 'Metafor', 'oo_metaB', 'bin'))
    addPath(os.path.join(topdir, 'Metafor', 'oo_meta'))
    addPath(os.path.join(topdir, 'Metafor', 'linuxbin'))
    addPath(os.path.join(topdir, 'NativeSolid', 'bin'))
    addPath(os.path.join(topdir, 'ModalSolver'))
    addPath(os.path.join(topdir, 'waves'))
    addPath(os.path.join(topdir, 'PFEM'))
    addPath(os.path.join(topdir, 'SU2', 'bin'))
    print 'PYTHONPATH = %s\n' % sys.path

def main():
    # Global variables
    global __file__

    # Find solvers and set paths
    setPath()

    # Parse options
    import os, argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('file', nargs=1, help='Python input file')
    parser.add_argument('-n', dest='n', help='Number of threads or processes', default=1, type=int)
    parser.add_argument('--nogui', action='store_true', help='Disable GUI', default=False)
    args = parser.parse_args()

    # Process
    for file in args.file:
        if not os.path.isfile(file):
            raise RuntimeError('Could not find file:', file)
        else:
            file = os.path.abspath(file)
            # change dir
            import cupydo.utilities as cutil
            cutil.setDirs(file)
            # split streams
            __file__ = file
            print "[run.py] __file__", __file__
            tee = Tee('stdout.txt')
            # run
            import time, platform
            print '-' * 80
            print 'Starting test:', file
            print 'Time:', time.strftime('%c')
            print 'Hostname:', platform.node()
            execfile(file, globals(), globals())

if __name__ == '__main__':
    main()