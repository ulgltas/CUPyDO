// Copyright 2018 University of Liege
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

// Gmsh project created on Wed Sep 14 10:38:00 2016
L = 0.03;
h = 0.001;
a = 0.0034;
b = a;
dx = 0.0117;
dy = 0.0243;
N = 10;

Point(1) = {0, 0, 0, 1.0};
Point(2) = {h, 0, 0, 1.0};
Point(3) = {h, L, 0, 1.0};
Point(4) = {0, L, 0, 1.0};
Point(5) = {-dx, dy, 0, 1.0};
Point(6) = {-dx+a, dy, 0, 1.0};
Point(7) = {-dx+a, dy+b, 0, 1.0};
Point(8) = {-dx, dy+b, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

Line Loop(9) = {1, 2, 3, 4};
Plane Surface(10) = {9};
Line Loop(11) = {5, 6, 7, 8};
Plane Surface(12) = {11};

Physical Line(13) = {1};          // Clamping
Physical Line(14) = {5, 6, 7, 8}; // Free surface
Physical Line(15) = {4};	  // FSI interface
Physical Surface(16) = {10};	  // Beam
Physical Surface(17) = {12};      // Bird 

Transfinite Line{1} = (h/a)*N+2;
Transfinite Line{2} = (3/2)*(L/a)*N+1;
Transfinite Line{3} = (h/a)*N+2;
Transfinite Line{4} = (3/2)*(L/a)*N+1;
Transfinite Line{5} = N+1;
Transfinite Line{6} = N+1;
Transfinite Line{7} = N+1;
Transfinite Line{8} = N+1;
Transfinite Surface{10} = {1,2,3,4};
Recombine Surface{10};
