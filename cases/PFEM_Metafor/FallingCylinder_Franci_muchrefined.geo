N = 200;
Ncurve = 48;
l = 0.02;
a = 0.0025;
h0 = 3.0*l;

Point(1) = {0, 0, 0, 1.0};
Point(2) = {2*l, 0, 0, 1.0};
Point(3) = {2*l, 30*a, 0, 1.0};
Point(4) = {0, 30*a, 0, 1.0};
Point(5) = {l, h0-a, 0, 1.0};
Point(6) = {l, h0+a, 0, 1.0};
Point(7) = {l, h0, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Circle(5) = {5, 7, 6};
Circle(6) = {6, 7, 5};

Line Loop(7) = {1, 2, 3, 4};
Line Loop(8) = {5, 6};
Plane Surface(9) = {7, 8};
Plane Surface(10) = {8};

Physical Line(11) = {1};
Physical Line(12) = {2};
Physical Line(13) = {3};
Physical Line(14) = {4};
Physical Line(15) = {6, 5};
Physical Surface(16) = {9};
Physical Surface(17) = {10};

Transfinite Line {1} = N;
Transfinite Line {2} = 3*N/2;
Transfinite Line {3} = N;
Transfinite Line {4} = 3*N/2;
Transfinite Line {5} = Ncurve;
Transfinite Line {6} = Ncurve;