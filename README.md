# CUPyDO
FSI tools for partinioned coupling between generic solid and fluid solvers.
[![Apache License Version 2.0](https://img.shields.io/badge/license-Apache_2.0-green.svg)](LICENSE)

## Features
As of December 2018, interfaces for the following solvers are implemented:

- Solid:
  - Metafor --> A Nonlinear Finite Element solid solver developed at the University of Liège (http://metafor.ltas.ulg.ac.be/dokuwiki/)
  - RBM --> A dynamic 2dof pitch/plunge solid solver developed at the University of Liège (https://github.com/ulgltas/NativeSolid)
  - Modal Solver -> A static/dynamic modal solver developed at the University of Liège (https://github.com/ulgltas/ModalSolver)
  - GetDP --> A free finite element software and a general environment for the treatment of discrete problems, developed at University of Liège (http://getdp.info/)
- Fluid:
  - PFEM --> Particle Finite Element Method fluid solver developed at the University of Liège (https://github.com/ulgltas/PFEM)
  - SU2 --> Open-source CFD code developed at Stanford University (http://su2.stanford.edu/)
  - Flow --> A Full Potential Finite Element fluid solver, part of the waves project, developed at the University of Liège (https://github.com/ulgltas/waves)

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
make install
```
Notes:
* MPI should be enabled if you want to use the MPI version of SU2. Keep this option disabled otherwise. 
* `path/to/petsc/include` is usually `/usr/lib/petscdir/version/` on Ubuntu/Debian
* The "install" step is mandatory. It copies the binaries in `CUPyDO/ccupydo`.


Run test battery
```bash
ctest [-R include_pattern] [-E exclude_pattern] -j4
```

Run test file
```bash
python run.py path/to/testfile_fsi.py -n4
```

## Interfaced solvers compilation (linux - gcc)
Brief instructions to compile interfaced solvers. The full documentation is available on the aforementioned websites.
The directories containing the external solvers must be placed at the same level as CUPyDO directory.

### Common packages
```bash
sudo apt-get install gmsh
sudo apt-get install libtbb-dev
sudo apt-get install libvtk6.3 libvtk6-dev libvtk6-qt-dev python-vtk6 python-pyqt5
```

### Metafor
[Linux, Windows, macOS]

Required packages
```bash
sudo apt install subversion
sudo apt install bison flex 
mkdir Metafor && cd Metafor
git clone git@github.com:ulgltas/linuxbin.git
```
Compilation
```bash
svn co svn+ssh://username@blueberry.ltas.ulg.ac.be/home/metafor/SVN/oo_meta/trunk oo_meta
mkdir oo_metaB && cd oo_metaB
cmake -C ../oo_meta/CMake/configMachine-CUPyDO.cmake [-DCMAKE_INSTALL_PREFIX=/path/to/Metafor/install/folder] [-DCMAKE_BUILD_TYPE=Debug] ../oo_meta
make -j4
```
Notes: 
* Metafor cannot be built with the "parasolid" interface which uses the same class names as `waves/fwk`
* A "sutdent" configuration can be used but some tests will not pass due to license restrictions.

### RBM [NativeSolid]
[Linux]

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
Windows: FIX LAPACKE


### SU2
[Linux]
```bash
sudo apt-get install autoconf
git clone git@github.com:su2code/SU2.git
cd SU2
git checkout tags/v6.2.0
unset MKLROOT   # <= MKL should be disabled
./bootstrap
./configure --prefix=/path/to/SU2/install/folder CXXFLAGS="-O3" [--enable-mpi --with-cc=/path/to/mpicc --with-cxx=/path/to/mpicxx] --enable-PY_WRAPPER [--enable-tecio]
make -j4
make install
```
Notes:
* The INSTALL step is mandatory!
* MPI should be enabled/disabled in both CUPyDO and SU2.


### GetDP
Todo

### Waves/flow
[Linux, Windows, macOS]

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
cmake -C disable-trilinos.cmake ..
make -j4
```
Note:
* `disable-trilinos.cmake` is not mandatory but it helps save much build time.

### PFEM
[Linux, Windows, macOS]
```bash
git clone git@github.com:ulgltas/waves.git 
[build waves as above]
git@github.com:ulgltas/PFEM.git
cd PFEM
mkdir build && cd build
cmake ..
make -j4
```


## Examples
Examples of simulations are available in [CUPyDO/tests](https://github.com/ulgltas/CUPyDO/tree/master/tests) and [CUPyDO/cases](https://github.com/ulgltas/CUPyDO/tree/master/cases).

![Screenshot](/tests/fsi_examples.png)

## Publications:
Cerquaglia M.L., Thomas D., Boman R., Terrapon V.E., Ponthot J.-P., [A fully partitioned Lagrangian framework for FSI problems characterized by free surfaces, large solid deformations and displacements, and strong added-mass effects](https://doi.org/10.1016/j.cma.2019.01.021), Computer Methods in Applied Mechanics and Engineering, in press (2019)

Thomas D., Cerquaglia M.L., Boman R., Economon T.D., Alonso J.J., Dimitriadis G., Terrapon V.E., [CUPyDO - An integrated Python environment for coupled multi-physics simulations](https://doi.org/10.1016/j.advengsoft.2018.05.007), Advances in Engineering Software 128:69-85 (2019)
