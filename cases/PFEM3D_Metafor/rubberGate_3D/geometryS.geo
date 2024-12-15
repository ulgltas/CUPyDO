H = 0.079;
S = 0.005;
L1 = 0.1;
W = 0.1;

eps = 1e-5;

P = 3;
M = 41;
N = 32;

// Point List

Point(1) = {L1, 0, H};
Point(2) = {L1+S, 0, H};
Point(3) = {L1, W, H};
Point(4) = {L1+S, W, H};
Point(5) = {L1, eps, eps};
Point(6) = {L1+S, eps, eps};
Point(7) = {L1, W-eps, eps};
Point(8) = {L1+S, W-eps, eps};

// Line List

Line(1) = {1, 2};
Line(2) = {3, 4};
Line(3) = {2, 4};
Line(4) = {1, 3};


Line(5) = {5, 6};
Line(6) = {7, 8};
Line(7) = {6, 8};
Line(8) = {5, 7};
Line(9) = {5, 1};
Line(10) = {6, 2};
Line(11) = {8, 4};
Line(12) = {7, 3};

Transfinite Line{1} = P;
Transfinite Line{2} = P;
Transfinite Line{5} = P;
Transfinite Line{6} = P;

Transfinite Line{3} = M;
Transfinite Line{4} = M;
Transfinite Line{7} = M;
Transfinite Line{8} = M;

Transfinite Line{9} = N;
Transfinite Line{10} = N;
Transfinite Line{11} = N;
Transfinite Line{12} = N;

// Mesh of the Solid

Curve Loop(1) = {1, 3, -2, -4};
Curve Loop(2) = {8, 6, -7, -5};
Curve Loop(3) = {5, 10, -1, -9};
Curve Loop(4) = {12, 2, -11, -6};
Curve Loop(5) = {7, 11, -3, -10};
Curve Loop(6) = {9, 4, -12, -8};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};

Transfinite Surface{1:6};
Recombine Surface{1:6};

Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

Transfinite Volume(1);
Recombine Volume(1);

// Physical Surfaces

Physical Volume("Solid") = {1};
Physical Surface("FSI") = {6};
Physical Surface("Side") = {3, 4};
Physical Surface("Base") = {1};

Mesh 3;