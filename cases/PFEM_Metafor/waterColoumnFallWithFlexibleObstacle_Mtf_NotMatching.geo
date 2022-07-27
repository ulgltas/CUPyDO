// Gmsh project created on Sat Aug 06 14:21:51 2016
L = 0.146;
a = 0.012;
b = (20.0/3.0)*a;
N = 10;

Point(5) = {2*L, 0, 0, 1.0};
Point(6) = {2*L, b, 0, 1.0};
Point(7) = {2*L+a, b, 0, 1.0};
Point(8) = {2*L+a, 0, 0, 1.0};

Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 8};
Line(11) = {8, 5};

Line Loop(14) = {6, 7, 8, 11};
Plane Surface(15) = {14};

Physical Line(17) = {6, 7, 8};        //Obstacle surface: FSI interface
Physical Line(19) = {11};             //Obstacle base
Physical Surface(21) = {15};          //Obstacle

Transfinite Line{6} = 2*(b/L)*N;
Transfinite Line{7} = (a/L)*N+3;
Transfinite Line{8} = 2*(b/L)*N;
Transfinite Line{11} = (a/L)*N+3;
Transfinite Surface{15} = {5,8,7,6};
Recombine Surface{15};