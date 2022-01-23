d = 0.1;
N = 80;

R = 2.25;
H = 3.75;
B = 4.87;
h = 2.5;
b = 1.3;
s = 0.2;

// Points List

Point(1) = {-B/2,H+h,0,d};
Point(2) = {-b/2,H,0,d};
Point(3) = {-R,H,0,d};
Point(4) = {-(R+s),H,0,d};
Point(5) = {-R,0,0,d};
Point(6) = {-(R+s),0,0,d};
Point(7) = {B/2,H+h,0,d};
Point(8) = {b/2,H,0,d};
Point(9) = {R,H,0,d};
Point(10) = {R+s,H,0,d};
Point(11) = {R,0,0,d};
Point(12) = {R+s,0,0,d};
Point(13) = {0,0,0,d};

// Lines List

Line(1) = {3,4};
Line(2) = {4,6};
Line(3) = {6,5};
Line(4) = {5,3};

Line(5) = {11,12};
Line(6) = {12,10};
Line(7) = {10,9};
Line(8) = {9,11};

Circle(9) = {6,13,12};
Circle(10) = {5,13,11};

Line(11) = {1,2};
Line(12) = {2,8};
Line(13) = {8,7};
Line(14) = {7,1};

Line(15) = {2,3};
Line(16) = {8,9};
Line(17) = {1,4};
Line(18) = {7,10};

// Solid Surface

Curve Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
Transfinite Surface{1};
Recombine Surface{1};

Curve Loop(2) = {5,6,7,8};
Plane Surface(2) = {2};
Transfinite Surface{2};
Recombine Surface{2};

Transfinite Line{9} = N;
Transfinite Line{10} = N;

Curve Loop(3) = {-3,9,-5,-10};
Plane Surface(3) = {3};
Transfinite Surface{3};
Recombine Surface{3};

Physical Surface("Solid") = {1,2,3};

// Fluid Surface

Curve Loop(4) = {11,12,13,14};
Plane Surface(4) = {4};
Physical Surface("Fluid") = {4};

// Physical Boundaries

Physical Curve("FSInterface") = {1,2,4,6,7,8,9,10};
Physical Curve("Reservoir") = {11,13,15,16};
Physical Curve("FreeSurface") = {12,14};
Physical Curve("SolidBase") = {1,7};

Physical Curve("PolyL") = {1,15,11,17};
Physical Curve("PolyR") = {7,16,13,18};

Mesh 2;