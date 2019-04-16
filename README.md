# CUPyDO
FSI tools for partinioned coupling between generic solid and fluid solvers.

## Features
As of December 2018, interfaces for the following solvers are implemented:

- Solid:
  - Metafor --> A Nonlinear Finite Element solid solver developed at the University of Liege (http://metafor.ltas.ulg.ac.be/dokuwiki/)
  - RBM --> A dynamic 2dof pitch/plunge solid solver developed at the University of Liege (https://github.com/ulgltas/NativeSolid)
- Fluid:
  - PFEM --> Particle Finite Element Method fluid solver developed at the University of Liege (http://metafor.ltas.ulg.ac.be/dokuwiki/poliflows/start)
  - SU2 --> Open-source CFD code developed at Stanford University (http://su2.stanford.edu/)
  - GetDP --> A free finite element software and a general environment for the treatment of discrete problems, developed at University of Liege (http://getdp.info/)
  - Flow --> A Full Potential Finite Element fluid solver, part of the waves project, developed at the University of Liege ()

## CUPyDO Compilation (linux - gcc)
Required packages
```bash
sudo apt-get install build-essential
sudo apt-get install cmake
sudo apt-get install python2.7 python2.7-dev libpython2.7-dev
sudo apt-get install swig
sudo apt-get install python-numpy python-scipy

```
Optional packages (parallel build, required for SU2+MPI)
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
[export INCLUDE=${INCLUDE}:/path/to/petsc/include]
cmake [-DWITH_MPI=ON] [-DCMAKE_BUILD_TYPE=Debug] ..
make -j4
```
*NB: `path/to/petsc/include` is usually `/usr/lib/petscdir/version/` on ubutun/debian*

Install
```bash
make install
make clean
```

Run test battery
```bash
./run.sh all 4
```

Run test file
```bash
./run.sh path/to/testfile 4
```

## Interfaced solvers compilation (linux - gcc)
Brief instructions to compile interfaced solvers. The full documentation is available on the aforementioned websites.
The directories containing the external solvers must be placed next to CUPyDO directory.

### Common packages
```bash
sudo apt-get install gmsh
sudo apt-get install libtbb-dev
sudo apt-get install libvtk6.3 libvtk6-dev libvtk6-qt-dev python-vtk6 python-pyqt5
```

### Metafor
Required packages
```bash
sudo apt install subversion
sudo apt install bison flex 
mkdir Metafor_Home && cd Metafor_Home
git clone git@github.com:ulgltas/linuxbin.git
```
Compilation
```bash
svn co svn+ssh://username@blueberry.ltas.ulg.ac.be/home/metafor/SVN/oo_meta/trunk oo_meta
mkdir oo_metaB && cd oo_metaB
cmake -C ../oo_meta/CMake/configMachine-CUPyDO.cmake <-DCMAKE_INSTALL_PREFIX=/path/to/Metafor/install/folder> <-DCMAKE_BUILD_TYPE=Debug> ../oo_meta
make -j4
<make install>
```

### RBM
Required packages
```bash
sudo apt-get install liblapacke-dev
sudo apt-get install libatlas-base-dev
```
Compilation
```bash
git clone git@github.com:ulgltas/NativeSolid.git
cd NativeSolid
mkdir build && cd build
cmake ..
make -j4
```

### PFEM
Todo

### SU2
```bash
sudo apt-get install autoconf
git clone git@github.com:su2code/SU2.git
cd SU2
git checkout tags/v6.2.0
./bootstrap
./configure --prefix=/path/to/SU2/install/folder CXXFLAGS="-O3" --enable-mpi --with-cc=/path/to/mpicc --with-cxx=/path/to/mpicxx --enable-PY_WRAPPER <--enable-tecio>
make -j4
make install
make clean
```

NOTE: Romain/garfield:
```
./configure --prefix=/home/boman/dev/CUPyDO/SU2 CXXFLAGS="-O3" --enable-mpi --with-cc=mpicc --with-cxx=mpicxx --enable-PY_WRAPPER --enable-tecio
```


### GetDP
Todo

### Waves/flow
Required packages
```bash
sudo apt-get install libgmm++-dev libeigen3-dev
sudo apt-get install libmumps-seq-dev libopenblas-dev
```
Optionnal packages
```bash
sudo apt-get install python-matplotlib
```
Compilation
```bash
git clone git@github.com:ulgltas/waves.git
cd waves
mkdir build && cd build
cmake -C ../CMake/ubuntu.cmake ..
make -j4
```

## Examples
Examples of simulations available in /tests.

![Screenshot](/tests/fsi_examples.png)

## Publications:
Cerquaglia M.L., Thomas D., Boman R., Terrapon V.E., Ponthot J.-P., [A fully partitioned Lagrangian framework for FSI problems characterized by free surfaces, large solid deformations and displacements, and strong added-mass effects](https://doi.org/10.1016/j.cma.2019.01.021), Computer Methods in Applied Mechanics and Engineering, in press (2019)

Thomas D., Cerquaglia M.L., Boman R., Economon T.D., Alonso J.J., Dimitriadis G., Terrapon V.E., [CUPyDO - An integrated Python environment for coupled multi-physics simulations](https://doi.org/10.1016/j.advengsoft.2018.05.007), Advances in Engineering Software 128:69-85 (2019)

[![Apache License Version 2.0](https://img.shields.io/badge/license-Apache_2.0-green.svg)](LICENSE)

