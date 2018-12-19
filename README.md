# CUPyDO
FSI tools for partinioned coupling between generic solid and fluid solvers.

## Features
As of December 2018, interfaces for the following solvers are implemented:

- Solid:
  - Metafor --> A Nonlinear Finite Element solid solver developed at the University of Liège (http://metafor.ltas.ulg.ac.be/dokuwiki/)
  - RBM --> A dynamic 2dof pitch/plunge solid solver developed at the University of Liège (https://github.com/ulgltas/NativeSolid)
- Fluid:
  - PFEM --> Particle Finite Element Method fluid solver developed at the University of Liège (http://metafor.ltas.ulg.ac.be/dokuwiki/poliflows/start)
  - SU2 --> Open-source CFD code developed at Stanford University (http://su2.stanford.edu/)
  - GetDP --> A free finite element software and a general environment for the treatment of discrete problems, developed at University of Liège (http://getdp.info/)
  - Flow --> A Full Potential Finite Element fluid solver, part of the waves project, developed at the University of Liège ()

## CUPyDO Compilation (linux - gcc)
Required packages
```bash
sudo apt-get install build-essential
sudo apt-get install cmake
sudo apt-get install python2.7 python2.7-dev libpython2.7-dev
sudo apt-get install swig
sudo apt-get install python-numpy python-scipy

```
Optional packages (parallel build, required for SU2)
```bash
sudo apt-get install libopenblas-dev liblapack-dev
sudo apt-get install openmpi-bin openmpi-common libopenmpi-dev
sudo apt-get install petsc-dev
sudo apt-get install python-mpi4py
sudo apt-get install python-petsc4py
```
Compilation
```bash
mkdir build && cd build
<INCLUDE=${INCLUDE}:/path/to/petsc/include>
cmake <-DWITH_MPI=ON> <-DCMAKE_BUILD_TYPE=Release> ..
make -j4
```
Install
```bash
make install
make clean
```
Note that the path to executables of the interfaced solvers should be put in your python path.

## Interfaced solvers compilation (linux - gcc)
Brief instructions to compile interfaced solvers. The full documentation is available on the aforementioned websites.

### Metafor
Todo

### RBM
Required packages
```bash
sudo apt-get install liblapacke-dev
sudo apt-get install liblatlas-base-dev
```
Compilation
```bash
git clone https://github.com/ulgltas/NativeSolid.git
cd NativeSolid
mkdir build && cd build
cmake ..
make -j4
```

### PFEM
Todo

### SU2
```bash
git clone https://github.com/su2code/SU2.git
cd SU2
git checkout tags/v6.0.0
./configure --prefix=/path/to/SU2/install/folder CXXFLAGS="-O3" --enable-mpi --with-cc=/path/to/mpicc --with-cxx=/path/to/mpicxx --enable-PY_WRAPPER <--enable-tecio>
make -j4
make install
make clean
```

### GetDP
Todo

### Flow
Required packages
```bash
sudo apt-get install libgmm++-dev libeigen3-dev
sudo apt-get install libtbb-dev
sudo apt-get install libmumps-seq-dev libopenblas-dev
sudo apt-get install gmsh
sudo apt-get install python-pyqt5
```
Optionnal packages
```bash
sudo apt-get install python-vtk6 libvtk6.3
sudo apt-get install python-matplotlib
```
Compilation
```bash
git clone
cd waves
mkdir build && cd build
cmake -C ../CMake/ubuntu.cmake ..
make -j4
```

### Environment variables
```bash
export PYTHONPATH="${PYTHONPATH}:/path/to/NativeSolid/bin"
export PYTHONPATH="${PYTHONPATH}:/path/to/SU2/install/folder/bin"
export PYTHONPATH="${PYTHONPATH}:/path/to/waves/"
export PYTHONPATH="${PYTHONPATH}:/path/to/waves/build/bin"
```

## Examples
Examples of simulations available in /tests.

![Screenshot](/tests/fsi_examples.png)

[![Apache License Version 2.0](https://img.shields.io/badge/license-Apache_2.0-green.svg)](LICENSE)

