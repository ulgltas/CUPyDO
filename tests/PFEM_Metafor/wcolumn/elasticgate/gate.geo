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

// Gmsh project created on Wed Aug 24 16:48:21 2016
A = 0.1;   //m
H = 0.14;  //m
L = 0.079; //m
s = 0.005; //m
Nel = 30;
N = Nel+1;
d = 0.; // A/(2*Nel);   //m

// Point(1) = {0, H, 0, 1.0};
Point(2) = {0, 0, 0, 1.0};
Point(3) = {5*A, 0, 0, 1.0};
Point(4) = {6*A, 0, 0, 1.0};
Point(5) = {6*A, H, 0, 1.0};
Point(6) = {5*A, H, 0, 1.0};
Point(7) = {5*A, L, 0, 1.0};
Point(8) = {5*A-s, L, 0, 1.0};
Point(9) = {5*A-s, d, 0, 1.0};
Point(10) = {5*A, d, 0, 1.0};

// Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 10};
Line(8) = {10, 3};
Line(9) = {7, 8};
Line(10) = {8, 9};
Line(11) = {9, 10};

Line Loop(12) = {3, 4, 5, 6, 7, 8};
Plane Surface(13) = {12};
Line Loop(14) = {9, 10, 11, -7};
Plane Surface(15) = {14};

Physical Line("Reservoir",16) = {3, 4, 6}; //Reservoir
Physical Line("FSInterface",17) = {7};        //Obstacle surface: FSI interface
Physical Line("FreeSurface",18) = {5, 8};              //Free surface
Physical Line("ObstacleBase",19) = {9};             //Obstacle base
Physical Surface("WaterColumn",20) = {13};          //Water column
Physical Surface("Obstacle",21) = {15};          //Obstacle
Physical Line("Outside",22) = {2};        //Outside part of the reservoir
Physical Surface("Contact",100) = {};          //Solid-solid contact

// Transfinite Line{1} = (H/A)*N;
Transfinite Line{2} = 5*N;
Transfinite Line{3} = N;
Transfinite Line{4} = (H/A)*N;
Transfinite Line{5} = N;
Transfinite Line{6} = ((H-L)/A)*N+1;
Transfinite Line{7} = ((L-d)/A)*N;
Transfinite Line{8} = (d/A)*N;
Transfinite Line{9} = (s/A)*N+3;
Transfinite Line{10} = ((L-d)/A)*N;
Transfinite Line{11} = (s/A)*N+3;
Transfinite Surface{15} = {10, 7, 8, 9};
Recombine Surface{15};