R = 0.06;
L = 0.4;
h = 0.00635;

N = 10;
nNodes = N+1;
d = 2.5*(R/N);

// Points List

Point(1) = {0,-h,0,1.0};
Point(2) = {L,-h,0,1.0};
Point(3) = {L,0,0,1.0};
Point(4) = {0,0,0,1.0};

// Lines List

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

// Fluid Surface

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

Transfinite Line{1} = 1.5*(L/R)*N+1;
Transfinite Line{2} = 2*(h/R)*N+1;
Transfinite Line{3} = 1.5*(L/R)*N+1;
Transfinite Line{4} = 2*(h/R)*N+1;

Transfinite Surface{1};
Recombine Surface{1};

// Physical Boundaries

Physical Surface("Solid",11) = {1};
Physical Line("FSInterface",12) = {3};
Physical Line("Clamped",13) = {2};
Physical Line("Axis",14) = {4};

Mesh 2;