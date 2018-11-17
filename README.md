# CUPyDO
FSI tools for partinioned coupling between generic solid and fluid solvers.

## Features
For the moment, interfaces for the following solvers are implemented:
    
- Metafor --> A Nonlinear Finite Element solid solver developed at the University of Liège (http://metafor.ltas.ulg.ac.be/dokuwiki/)
- PFEM --> Particle Finite Element Method fluid solver developed at the University of Liège (http://metafor.ltas.ulg.ac.be/dokuwiki/poliflows/start)
- SU2 --> Open-source CFD code developed at Stanford University (http://su2.stanford.edu/)
- GetDP --> A free finite element software and a general environment for the treatment of discrete problems, developed at University of Liège (http://getdp.info/)
- Flow --> A Full Potential Finite Element fluid solver developed at the University of Liège

Examples of simulations available in /tests.

## Compilation (linux - gcc)
Required packages
```bash
sudo apt-get install cmake
sudo apt-get install python-dev
sudo apt-get install swig
sudo apt-get install python-numpy
```
Optional packages
```
python-mpi4py
```
Compilation
```
mkdir build && cd build
make -j4
```
Install
```
make install
```

![Screenshot](/tests/fsi_examples.png)

[![Apache License Version 2.0](https://img.shields.io/badge/license-Apache_2.0-green.svg)](LICENSE)

