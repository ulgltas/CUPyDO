L = 0.146;
s = 0.012;
h = 0.08;
d = 0.015;

N = 10;
M = 3;

// Points List

Point(1) = {2*L,0,0,d};
Point(2) = {2*L,h,0,d};
Point(3) = {2*L+s,h,0,d};
Point(4) = {2*L+s,0,0,d};

// Lines List

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {1,4};

// Solid Surface

Transfinite Line{1} = N;
Transfinite Line{3} = N;
Transfinite Line{2} = M;
Transfinite Line{4} = M;

Curve Loop(1) = {-3,-2,-1,4};
Plane Surface(1) = {1};

Transfinite Surface{1};
Recombine Surface{1};

// Physical Boundaries

Physical Curve("FSI") = {1,2,3};
Physical Curve("SolidBase") = {4};
Physical Surface("Solid") = {1};

Mesh 2;