R = 0.06;
L = 0.4;
h = 0.00635;

N = 10;
nNodes = N+1;
D = 2.5*(R/N);


d = 0.0054;

// Points List

Point(1) = {0,D,0,d};
Point(2) = {0,D+R,0,d};
Point(3) = {R,D+R,0,d};
Point(4) = {R,D+R+2*R,0,d};
Point(5) = {0,D+R+2*R,0,d};
Point(6) = {0,0,0,d};
Point(7) = {L,0,0,d};

// Lines List

Circle(1) = {1,2,3};
Line(2) = {3,4};
Line(3) = {4,5};
Line(4) = {5,1};
Line(5) = {6,7};

// Fluid Surface

Line Loop(1) = {1,2,3,4};
Line Loop(2) = {1,2,3,4};

Plane Surface(1) = {1};
Plane Surface(2) = {2};

// Transfinite Line{1} = 1.5*(N+1);
// Transfinite Line{2} = 2*(N+1);
// Transfinite Line{3} = (N+1);
// Transfinite Line{4} = 3*(N+1);
Transfinite Line{5} = 1.5*(L/R)*N+1;

// Physical Boundaries

Physical Surface("Fluid") = {1};
Physical Line("FSInterface") = {5};
Physical Line("FreeSurface") = {1,2,3};

Mesh 2;