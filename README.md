# Copyright 2018 University of Liège
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# CUPyDO
FSI tools for partinioned coupling between generic solid and fluid solvers.

For the moment, interfaces for the following solvers are implemented:
    
- Metafor --> A Nonlinear Finite Element solid solver developed at the University of Liège (http://metafor.ltas.ulg.ac.be/dokuwiki/)
- PFEM --> Particle Finite Element Method fluid solver developed at the University of Liège (http://metafor.ltas.ulg.ac.be/dokuwiki/poliflows/start)
- SU2 --> Open-source CFD code developed at Stanford University (http://su2.stanford.edu/)
- GetDP --> A free finite element software and a general environment for the treatment of discrete problems, developed at University of Liège (http://getdp.info/)

Examples of simulations available in /tests.

![Screenshot](/tests/fsi_examples.png)

