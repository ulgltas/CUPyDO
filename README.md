# CUPyDO
FSI tools for partinioned coupling between generic solid and fluid solvers.

[![Apache License Version 2.0](https://img.shields.io/badge/license-Apache_2.0-green.svg)](LICENSE)

## Solvers
CUPyDO currently features interfaces for the following solvers:
- **Solid**
  - Metafor **[v3491]** (http://metafor.ltas.ulg.ac.be/dokuwiki/)
  - A Nonlinear Finite Element solid solver developed at the University of Liège.
  ---
  - RBM **[X]** (https://github.com/ulgltas/NativeSolid)
  - A dynamic 2dof pitch/plunge solid solver developed at the University of Liège.
  ---
  - SU2 **[X]** (https://su2code.github.io/)
  - Open-source CFD code developed at Stanford University.
  ---
  - Modali **[v2.0]** (https://github.com/ulgltas/modali)
  - A static/dynamic modal solver developed at the University of Liège.
  ---
  - GetDP **[X]** (http://getdp.info/)
  - A free finite element software and a general environment for the treatment of discrete problems, developed at University of Liège.
  ---
  - pyBeam **[X]** (https://github.com/pyBeam/pyBeam)
  - A nonlinear beam finite element solver developed for aeronautical design applications.
  ---
- **Fluid**
  - PFEM **[v1.19]** (https://gitlab.uliege.be/am-dept/PFEM)
  - Particle Finite Element Method fluid solver developed at the University of Liège.
  ---
  - SU2 **[7.11\*]** (https://github.com/ulgltas/SU2/tree/fix_wrap_strong)
  - Open-source CFD code developed at Stanford University.
  ---
  - Flow **[v2.2.1]** (https://gitlab.uliege.be/am-dept/waves)
  - A Full Potential Finite Element fluid solver, part of the waves project, developed at the University of Liège.
  ---
  - VLM **[v2.0]** (https://github.com/ulgltas/VLM)
  - A Vortex Lattice Method, developed at the University of Liège.
  ---
  - PFEM3D **[Last]** (https://github.com/ImperatorS79/PFEM3D)
  - A 3D Particle Finite Element Method fluid solver developed at the University of Liège.

##  Features

Furthermore, CUPyDO features two interpolation alogrithms:
- Radial Basis Functions (RBF)
- Thin Plate Spline (TPS)

Finally, CUPyDO features two main couplers:
- Block-Gauss-Seidel (BGS) with
  - constant (static) relaxation
  - Aitken relaxation
- Interface Quasi-Newton with Inverse Least-Square (IQN_ILS)

## Compilation
Detailed build instructions can be found in the [wiki](https://github.com/ulgltas/CUPyDO/wiki/Installation).

## Examples
Examples of simulations are available in [CUPyDO/tests](https://github.com/ulgltas/CUPyDO/tree/master/tests) and [CUPyDO/cases](https://github.com/ulgltas/CUPyDO/tree/master/cases).

![Screenshot](/tests/fsi_examples.png)

## Publications:
Cerquaglia M.L., Thomas D., Boman R., Terrapon V.E., Ponthot J.-P., [A fully partitioned Lagrangian framework for FSI problems characterized by free surfaces, large solid deformations and displacements, and strong added-mass effects](https://doi.org/10.1016/j.cma.2019.01.021), Computer Methods in Applied Mechanics and Engineering, in press (2019)

Thomas D., Cerquaglia M.L., Boman R., Economon T.D., Alonso J.J., Dimitriadis G., Terrapon V.E., [CUPyDO - An integrated Python environment for coupled multi-physics simulations](https://doi.org/10.1016/j.advengsoft.2018.05.007), Advances in Engineering Software 128:69-85 (2019)
