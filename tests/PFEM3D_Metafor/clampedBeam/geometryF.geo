// Parameters

d = 0.03;
N = 31;

L = 1;
HF = 0.2;
HS = 0.02;

// Points List

Point(1) = {L,HS,0,d};
Point(2) = {0,HS,0,d};
Point(3) = {L,HS+HF,0,d};
Point(4) = {0,HS+HF,0,d};

// Lines List

Line(1) = {1,2};
Line(2) = {1,3};
Line(3) = {3,4};
Line(4) = {4,2};

// Fluid Surface

Transfinite Line{1} = N;
Transfinite Surface{1};
Recombine Surface{1};

Curve Loop(1) = {-1,2,3,4};
Plane Surface(1) = {1};
Physical Surface("Fluid") = {1};

// Boundaries

Physical Curve("FSInterface") = {1};
Physical Curve("FreeSurface") = {3};
Physical Curve("Wall") = {2,4};

Mesh 2;