// Copyright 2018 University of Li�ge
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
//     http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License. 

N = 10;
nNodes = N+1;
R = 0.06;
L = 0.4;
h = 0.00635;
d = 2.5*(R/N); 

Point(1) = {0, -h, 0, 1.0};
Point(2) = {L, -h, 0, 1.0};
Point(3) = {L, 0, 0, 1.0};
Point(4) = {0, 0, 0, 1.0};
Point(5) = {0, d, 0, 1.0};
Point(6) = {0, d+R, 0, 1.0};
Point(7) = {R, d+R, 0, 1.0};
Point(8) = {R, d+R+2*R, 0, 1.0};
Point(9) = {0, d+R+2*R, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Circle(5) = {5, 6, 7};
Line(6) = {7, 8};
Line(7) = {8, 9};
Line(8) = {9, 5};

Line Loop(9) = {5, 6, 7, 8};
Plane Surface(10) = {9};
Line Loop(11) = {1, 2, 3, 4};
Plane Surface(12) = {11};

Physical Line("FSInterface",13) = {3};          //FSI surface
Physical Line("Clamping",14) = {2, 4};       //Clamping + Axis of symmetry for the panel
Physical Line("Axis",15) = {8};	  //Axis of symmetry for the bird
Physical Line("FreeSurface",16) = {5, 6, 7}; 	  //Free surface
Physical Surface("Bird",17) = {10};      //Bird
Physical Surface("Panel",18) = {12};      //Panel

Transfinite Line{1} = 1.5*(L/R)*N+1;
Transfinite Line{2} = 2*(h/R)*N+1;
Transfinite Line{3} = 1.5*(L/R)*N+1;
Transfinite Line{4} = 2*(h/R)*N+1;
Transfinite Line{5} = 1.5*nNodes;
Transfinite Line{6} = 2*nNodes;
Transfinite Line{7} = nNodes;
Transfinite Line{8} = 3*nNodes;
Transfinite Surface{12} = {1,2,3,4};
Recombine Surface{12};