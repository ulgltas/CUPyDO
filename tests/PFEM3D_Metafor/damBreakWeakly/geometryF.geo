L = 0.146;
s = 0.012;
h = 0.08;
d = 0.015;

N = 10;
M = 3;

// Points List

Point(1) = {0,0,0,d};
Point(2) = {4*L+s,0,0,d};
Point(3) = {4*L,6*L,0,d};
Point(4) = {0,6*L,0,d};
Point(5) = {L,0,0,d};
Point(6) = {L,2*L,0,d};
Point(7) = {0,2*L,0,d};
Point(8) = {2*L,0,0,d};
Point(9) = {2*L,h,0,d};
Point(10) = {2*L+s,h,0,d};
Point(11) = {2*L+s,0,0,d};

// Lines List

Line(1) = {4,7};
Line(2) = {7,1};
Line(3) = {1,5};
Line(4) = {5,6};
Line(5) = {6,7};
Line(6) = {5,8};
Line(7) = {8,9};
Line(8) = {9,10};
Line(9) = {10,11};
Line(10) = {11,2};
Line(11) = {2,3};
Line(12) = {8,11};

// Fluid Surface

Transfinite Line{7} = N;
Transfinite Line{9} = N;
Transfinite Line{8} = M;
Transfinite Line{12} = M;

Curve Loop(1) = {2,3,4,5};
Plane Surface(1) = {1};

// Physical Boundaries

Physical Curve("Polytope") = {7,8,9,12};
Physical Curve("FSInterface") = {7,8,9};
Physical Curve("Reservoir") = {1,2,3,6,10,11};
Physical Curve("FreeSurface") = {5,4};
Physical Surface("Fluid") = {1};

Mesh 2;