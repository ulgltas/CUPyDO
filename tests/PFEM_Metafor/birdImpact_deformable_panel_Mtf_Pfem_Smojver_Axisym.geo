// Copyright 2018 University of Liège
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

N = 15;
nNodes = N+1;
R = 0.06;
L = 0.2032;
s = 0.0538;
h = 0.00635;
d = 2.5*(R/N); 

Point(1) = {0, -h, 0, 1.0};
Point(2) = {L, -h, 0, 1.0};
Point(3) = {L+s, -h, 0, 1.0};
Point(4) = {L+s, 0, 0, 1.0};
Point(5) = {L, 0, 0, 1.0};
Point(6) = {0, 0, 0, 1.0};
Point(7) = {0, d, 0, 1.0};
Point(8) = {0, d+R, 0, 1.0};
Point(9) = {R, d+R, 0, 1.0};
Point(10) = {R, d+R+2*R, 0, 1.0};
Point(11) = {0, d+R+2*R, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
Circle(7) = {7, 8, 9};
Line(8) = {9, 10};
Line(9) = {10, 11};
Line(10) = {11, 7};

Line Loop(11) = {7, 8, 9, 10};
Plane Surface(12) = {11}; //Bird
Line Loop(13) = {1, 2, 3, 4, 5, 6};
Plane Surface(14) = {13}; //Panel

Physical Line(15) = {2};       	  //Vertical fixation
Physical Line(16) = {3};       	  //Clamping
Physical Line(17) = {4, 5};       //FSI interface
Physical Line(18) = {6};       	  //Axis of symmetry for the panel
Physical Line(19) = {7, 8, 9}; 	  //Free surface of the bird
Physical Line(20) = {10};	  //Axis of symmetry for the bird
Physical Surface(21) = {12};      //Bird
Physical Surface(22) = {14};      //Panel

Transfinite Line{1} = 1.5*(L/R)*N+1;
Transfinite Line{2} = 1.5*(s/R)*N+1;
Transfinite Line{3} = 3*(h/R)*N+1;
Transfinite Line{4} = 1.5*(s/R)*N+1;
Transfinite Line{5} = 1.5*(L/R)*N+1;
Transfinite Line{6} = 3*(h/R)*N+1;
Transfinite Line{7} = 1.5*nNodes;
Transfinite Line{8} = 2*nNodes;
Transfinite Line{9} = nNodes;
Transfinite Line{10} = 3*nNodes;
Transfinite Surface{14} = {1,3,4,6};
Recombine Surface{14};