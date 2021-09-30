L = 0.15;
s = 0.015;
h = 0.07;
d = 0.01;

Point(1) = {0,0,0,d};
Point(2) = {4*L,0,0,d};
Point(3) = {4*L,4*L,0,d};
Point(4) = {0,3*L,0,d};
Point(5) = {L,0,0,d};
Point(6) = {L,2*L,0,d};
Point(7) = {0,2*L,0,d};
Point(8) = {2*L,0,0,d};
Point(9) = {2*L,h,0,d};
Point(10) = {2*L+s,h,0,d};
Point(11) = {2*L+s,0,0,d};

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

Curve Loop(1) = {2,3,4,5};
Plane Surface(1) = {1};
Transfinite Surface{1};

Curve Loop(2) = {-9,-8,-7,12};
Plane Surface(2) = {2};
Transfinite Surface{2};
Recombine Surface{2};

Physical Curve("FSInterface",16) = {7,8,9};
Physical Curve("Reservoir",17) = {1,2,3,6,10,11};
Physical Curve("FreeSurface",18) = {5,4};
Physical Curve("SolidBase",19) = {12};

Physical Surface("Fluid",20) = {1};
Physical Surface("Solid",21) = {2};

Mesh 2;