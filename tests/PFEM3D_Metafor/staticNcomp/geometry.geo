// Parameters

d = 0.03;
L = 1;
HF = 0.2;
HS = 0.02;

// Points List

Point(1) = {0,0,0,d};
Point(2) = {L,0,0,d};
Point(3) = {L,HS,0,d};
Point(4) = {0,HS,0,d};
Point(5) = {L,HS+HF,0,d};
Point(6) = {0,HS+HF,0,d};

// Lines List

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {3,5};
Line(6) = {5,6};
Line(7) = {6,4};

// Solid Surface

Curve Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
Physical Surface("Solid") = {1};
Transfinite Surface{1};
Recombine Surface{1};

// Fluid Surface

Curve Loop(2) = {-3,5,6,7};
Plane Surface(2) = {2};
Physical Surface("Fluid") = {2};

// Boundaries

Physical Curve("FSInterface") = {1,2,3,4};
Physical Curve("FreeSurface") = {6};
Physical Curve("Clamped") = {2,4};
Physical Curve("Wall") = {5,7};

// Builds 2D Mesh

Mesh 2;