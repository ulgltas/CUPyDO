# CUPyDO
FSI tools for partinioned coupling between generic solid and fluid solvers.

[![Apache License Version 2.0](https://img.shields.io/badge/license-Apache_2.0-green.svg)](LICENSE)

## Features
CUPyDO currently features interfaces for the following solvers:
- Solid:
  - Metafor --> A Nonlinear Finite Element solid solver developed at the University of Liège (http://metafor.ltas.ulg.ac.be/dokuwiki/)
  - RBM --> A dynamic 2dof pitch/plunge solid solver developed at the University of Liège (https://github.com/ulgltas/NativeSolid)
  - modali -> A static/dynamic modal solver developed at the University of Liège (https://github.com/ulgltas/modali)
  - GetDP --> A free finite element software and a general environment for the treatment of discrete problems, developed at University of Liège (http://getdp.info/)
- Fluid:
  - PFEM --> Particle Finite Element Method fluid solver developed at the University of Liège (https://gitlab.uliege.be/am-dept/PFEM)
  - SU2 --> Open-source CFD code developed at Stanford University (https://su2code.github.io/)
  - DART --> A Full Potential Finite Element fluid solver, part of the waves project, developed at the University of Liège (https://gitlab.uliege.be/am-dept/dartflo)
  - VLM --> A Vortex Lattice Method, developed at the University of Liège (https://github.com/ulgltas/VLM)

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
