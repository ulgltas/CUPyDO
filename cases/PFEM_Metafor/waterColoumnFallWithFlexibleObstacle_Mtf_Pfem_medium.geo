// Gmsh project created on Sat Aug 06 14:21:51 2016
L = 0.146;
a = 0.012;
b = (20.0/3.0)*a;
N = 45;

Point(1) = {0, 2*L, 0, 1.0};
Point(2) = {0, 0, 0, 1.0};
Point(3) = {L, 0, 0, 1.0};
Point(4) = {L, 2*L, 0, 1.0};
Point(5) = {2*L, 0, 0, 1.0};
Point(6) = {2*L, b, 0, 1.0};
Point(7) = {2*L+a, b, 0, 1.0};
Point(8) = {2*L+a, 0, 0, 1.0};
Point(9) = {4*L+a, 0, 0, 1.0};
Point(10) = {4*L+a, 6*L, 0, 1.0};
Point(11) = {0, 6*L, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {3, 5};
Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 8};
Line(9) = {8, 9};
Line(10) = {9, 10};
Line(11) = {8, 5};
Line(12) = {11, 1};

Line Loop(12) = {1, 2, 3, 4};
Plane Surface(13) = {12};
Line Loop(14) = {6, 7, 8, 11};
Plane Surface(15) = {14};

Physical Line(16) = {1, 2, 5, 9, 10, 12}; //Reservoir
Physical Line(17) = {6, 7, 8};        //Obstacle surface: FSI interface
Physical Line(18) = {3, 4};           //Free surface
Physical Line(19) = {11};             //Obstacle base
Physical Surface(20) = {13};          //Water column
Physical Surface(21) = {15};          //Obstacle
Physical Surface(100) = {};          //Solid-solid contact

Transfinite Line{1} = 2*N;
Transfinite Line{2} = N;
Transfinite Line{3} = 2*N;
Transfinite Line{4} = N;
Transfinite Line{5} = N;
Transfinite Line{9} = 3*N;
Transfinite Line{10} = (3/2)*6*N;
Transfinite Line{6} = 2*(b/L)*N;
Transfinite Line{7} = (a/L)*N+6;
Transfinite Line{8} = 2*(b/L)*N;
Transfinite Line{11} = (a/L)*N+6;
Transfinite Line{12} = 4*N;
Transfinite Surface{15} = {5,8,7,6};
Recombine Surface{15};