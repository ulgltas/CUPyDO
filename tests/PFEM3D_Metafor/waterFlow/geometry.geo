L = 0.1;
H2 = 0.1;
R = 0.01;
H1 = 0.05;
d = 0.006;

Point(1) = {-L,0,0,d};
Point(2) = {-2*R,0,0,d};
Point(3) = {-R,R,0,d};
Point(4) = {-R,H1,0,d};
Point(5) = {0,H1,0,d};
Point(6) = {R,H1,0,d};
Point(7) = {R,R,0,d};
Point(8) = {2*R,0,0,d};
Point(9) = {L,0,0,d};
Point(10) = {L,H1+H2,0,d};
Point(11) = {-L,H1+H2,0,d};
Point(12) = {-2*R,R,0,d};
Point(13) = {2*R,R,0,d};

Line(1) = {1,2};
Circle(2) = {2,12,3};
Line(3) = {3,4};
Circle(4) = {6,5,4};
Line(5) = {6,7};
Circle(6) = {7,13,8};
Line(7) = {8,9};
Line(8) = {9,10};
Line(9) = {10,11};
Line(10) = {11,1};
Line(11) = {2,8};

Curve Loop(1) = {1,2,3,-4,5,6,7,8,9,10};
Plane Surface(1) = {1};
Physical Surface("Fluid",100) = {1};

Curve Loop(2) = {11,-6,-5,4,-3,-2};
Plane Surface(2) = {2};
Physical Surface("Solid",200) = {2};

Physical Curve("FSInterface",300) = {2,3,4,5,6};
Physical Curve("Reservoir",400) = {1,7,9};
Physical Curve("SolidBase",600) = {11};
Physical Curve("Inlet",700) = {10};
Physical Curve("FreeSurface",800) = {8};

Mesh 2;