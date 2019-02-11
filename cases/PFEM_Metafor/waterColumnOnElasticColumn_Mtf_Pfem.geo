// Gmsh project created on Wed Jun 21 14:50:32 2017
Nel = 2;
N = Nel + 1;

Point(1) = {0, 0, 0, 1.0};
Point(2) = {0.05, 0, 0, 1.0};
Point(3) = {0.05, 0.25, 0, 1.0};
Point(4) = {0.05, 1, 0, 1.0};
Point(5) = {0, 1, 0, 1.0};
Point(6) = {0, 0.25, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 3};
Line(7) = {6, 1};

Line Loop(8) = {1, 2, -6, 7};
Plane Surface(9) = {8};
Line Loop(10) = {5, 6, 3, 4};
Plane Surface(11) = {10};

Physical Line(12) = {1};    // Clamped bottom
Physical Line(13) = {2, 7}; // Elastic column sides
Physical Line(14) = {6};    // FSI interface
Physical Line(15) = {3, 5}; // Water column sides
Physical Line(16) = {4};    // Top free surface

Physical Surface(17) = {9};  // Elastic column
Physical Surface(18) = {11}; // Water column

Transfinite Surface{9} = {1, 2, 3, 6};
Transfinite Surface{11} = {6, 3, 4, 5};
Transfinite Line{1} = N;
Transfinite Line{2} = 5*Nel+1; // 5*N+1;
Transfinite Line{3} = 15*Nel+1; // 15*N+1;
Transfinite Line{4} = N;
Transfinite Line{5} = 15*Nel+1; // 15*N+1;
Transfinite Line{6} = N;
Transfinite Line{7} = 5*Nel+1; // 5*N+1;