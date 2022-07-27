N = 100;

h2 = 2.5; 
h1 = 3.75;
R = 2.25;
b2 = 1.3;
b1 = 4.8714;
s = 0.2;

Point(1) = {-b1/2, h2, 0, 1.0};
Point(2) = {b1/2, h2, 0, 1.0};
Point(3) = {-b2/2, 0, 0, 1.0};
Point(4) = {b2/2, 0, 0, 1.0};
Point(5) = {-(R+s), 0, 0, 1.0};
Point(6) = {-R, 0, 0, 1.0};
Point(7) = {R, 0, 0, 1.0};
Point(8) = {R+s, 0, 0, 1.0};
Point(9) = {-(R+s), -h1, 0, 1.0};
Point(10) = {-R, -h1, 0, 1.0};
Point(11) = {0, -h1, 0, 1.0};
Point(12) = {R, -h1, 0, 1.0};
Point(13) = {R+s, -h1, 0, 1.0};

Line(1) = {2, 1};
Line(2) = {1, 3};
Line(3) = {3, 4};
Line(4) = {4, 2};
Line(5) = {5, 6};
Line(6) = {6, 3};
Line(7) = {4, 7};
Line(8) = {7, 8};
Line(9) = {6, 10};
Circle(10) = {10, 11, 12};
Line(11) = {12, 7};
Line(12) = {5, 9};
Circle(13) = {9, 11, 13};
Line(14) = {13, 8};
Line(15) = {9, 10};
Line(16) = {12, 13};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Line Loop(2) = {-5, 12, 15, -9};
Plane Surface(2) = {2};
Line Loop(3) = {-15, 13, -16, -10};
Plane Surface(3) = {3};
Line Loop(4) = {-8, -11, 16, 14};
Plane Surface(4) = {4};

Physical Line(1) = {1, 3}; 	// Top and bottom of the fluid
Physical Line(2) = {2, 4}; 	// Oblique walls
Physical Line(3) = {5, 8}; 	// Clamping of the elastic container
Physical Line(4) = {6, 7}; 	// Horizontal walls
Physical Line(5) = {9, 10, 11}; // Inner side of the elastic container: FSI interface

Physical Surface(6) = {1};      // Fluid
Physical Surface(7) = {2, 3, 4};// Elastic container
Physical Surface(100) = {};     // Solid-solid contact

Transfinite Line{1} = N;
Transfinite Line{2} = ((Sqrt(h2^2+(R-b2/2)^2))/b1)*N;
Transfinite Line{3} = (b2/b1)*N;
Transfinite Line{4} = ((Sqrt(h2^2+(R-b2/2)^2))/b1)*N;
Transfinite Line{5} = (s/b1)*N+1;
Transfinite Line{6} = ((R-b2/2)/b1)*N;
Transfinite Line{7} = ((R-b2/2)/b1)*N;
Transfinite Line{8} = (s/b1)*N+1;
Transfinite Line{9} = (h1/b1)*N;
Transfinite Line{10} = ((Pi*R)/b1)*N+1;
Transfinite Line{11} = (h1/b1)*N;
Transfinite Line{12} = (h1/b1)*N;
Transfinite Line{13} = ((Pi*R)/b1)*N+1;
Transfinite Line{14} = (h1/b1)*N;
Transfinite Line{15} = (s/b1)*N+1;
Transfinite Line{16} = (s/b1)*N+1;
Transfinite Surface{2} = {6, 5, 9, 10};
Transfinite Surface{3} = {10, 9, 13, 12};
Transfinite Surface{4} = {8, 7, 12, 13};
Recombine Surface{2};
Recombine Surface{3};
Recombine Surface{4};