// Gmsh project created on Wed Aug 24 16:48:21 2016
A = 0.1;   //m
H = 0.14;  //m
L = 0.079; //m
s = 0.005; //m
N = 30;
d = A/N;   //m

Point(1) = {0, H+d, 0, 1.0};
Point(2) = {0, 0, 0, 1.0};
Point(3) = {5*A, 0, 0, 1.0};
Point(4) = {6*A, 0, 0, 1.0};
Point(5) = {6*A, H+d, 0, 1.0};
Point(6) = {5*A, H+d, 0, 1.0};
Point(7) = {5*A, L+d, 0, 1.0};
Point(8) = {5*A-s, L+d, 0, 1.0};
Point(9) = {5*A-s, d, 0, 1.0};
Point(10) = {5*A, d, 0, 1.0};

Line(1) = {1, 2};
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

Physical Line(16) = {1, 2, 3, 4, 6}; //Reservoir
Physical Line(17) = {7};        //Obstacle surface: FSI interface
Physical Line(18) = {5, 8};              //Free surface
Physical Line(19) = {9};             //Obstacle base
Physical Surface(20) = {13};          //Water column
Physical Surface(21) = {15};          //Obstacle
Physical Surface(100) = {};          //Solid-solid contact

Transfinite Line{1} = (H/A)*N;
Transfinite Line{2} = 5*N;
Transfinite Line{3} = N;
Transfinite Line{4} = (H/A)*N;
Transfinite Line{5} = N;
Transfinite Line{6} = ((H-L)/A)*N+1;
Transfinite Line{7} = ((L-d)/A)*N;
Transfinite Line{8} = (d/A)*N;
Transfinite Line{9} = (s/A)*N+2;
Transfinite Line{10} = ((L-d)/A)*N;
Transfinite Line{11} = (s/A)*N+2;
Transfinite Surface{15} = {10, 7, 8, 9};
Recombine Surface{15};