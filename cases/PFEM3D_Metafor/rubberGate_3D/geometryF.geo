H = 0.079;
D = 0.14;
L1 = 0.1;
L2 = 0.1;
W = 0.1;

eps = 1e-5;
d = 3e-3;

M = 41;
N = 32;

// Point List

Point(1) = {0, 0, 0, d};
Point(2) = {L1, 0, 0, d};
Point(3) = {L1+L2, 0, 0, d};
Point(4) = {0, W, 0, d};
Point(5) = {L1, W, 0, d};
Point(6) = {L1+L2, W, 0, d};
Point(7) = {0, 0, D, d};
Point(8) = {L1, 0, D, d};
Point(9) = {L1+L2, 0, D, d};
Point(10) = {0, W, D, d};
Point(11) = {L1, W, D, d};
Point(12) = {L1+L2, W, D, d};
Point(13) = {L1, 0, H, d};
Point(14) = {L1, W, H, d};
Point(15) = {L1, eps, eps, d};
Point(16) = {L1, W-eps, eps, d};

// Line List

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {4, 5};
Line(4) = {5, 6};
Line(5) = {1, 4};
Line(6) = {3, 6};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {7, 10};
Line(14) = {1, 7};
Line(15) = {4, 10};
Line(16) = {3, 9};
Line(17) = {6, 12};
Line(18) = {13, 8};
Line(19) = {14, 11};
Line(20) = {13, 14};
Line(21) = {8, 11};
Line(22) = {2, 5};
Line(23) = {2, 13};
Line(24) = {5, 14};
Line(25) = {15, 16};
Line(26) = {15, 13};
Line(27) = {16, 14};
Line(28) = {2, 15};
Line(29) = {5, 16};

Transfinite Line{20} = M;
Transfinite Line{25} = M;

Transfinite Line{26} = N;
Transfinite Line{27} = N;

// Boundary Surface

Curve Loop(1) = {26, 20, -27, -25};
Curve Loop(2) = {18, 21, -19, -20};
Curve Loop(3) = {4, -6, -2, 22};
Curve Loop(4) = {4, 17, -11, -19, -24};
Curve Loop(5) = {2, 16, -8, -18, -23};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};

// Mesh of the Fluid

Curve Loop(6) = {7, 21, -10, -12};
Curve Loop(7) = {5, 3, -22, -1};
Curve Loop(8) = {1, 23, 18, -7, -14};
Curve Loop(9) = {15, 10, -19, -24, -3};
Curve Loop(10) = {5, 15, -12, -14};

Curve Loop(11) = {28, 26, -23};
Curve Loop(12) = {24, -27, -29};
Curve Loop(13) = {29, -25, -28, 22};

Plane Surface(6) = {6};
Plane Surface(7) = {7};
Plane Surface(8) = {8};
Plane Surface(9) = {9};
Plane Surface(10) = {10};

Plane Surface(11) = {11};
Plane Surface(12) = {12};
Plane Surface(13) = {13};

Surface Loop(1) = {1, 2, 6, 7, 8, 9, 10, 11, 12, 13};

Transfinite Surface{1};
Volume(1) = {1};

// Physical Entities

Physical Volume("Fluid") = {1};
Physical Surface("Reservoir") = {7, 8, 9, 10};
Physical Surface("Refine") = {2, 3, 4, 5};
Physical Surface("FSI") = {1};

// Mesh Size Field

Field[1] = Distance;
Field[1].SurfacesList = {1, 2, 3, 4, 5, 6};
Field[1].Sampling = 1e3;

Field[2] = MathEval;
Field[2].F = Sprintf("F1*0.1 + %g", d);

Field[3] = Box;
Field[3].VIn = 1e-2;
Field[3].VOut = 1e-2;
Field[3].XMin = 0;
Field[3].XMax = 1;
Field[3].YMin = 0;
Field[3].YMax = 1;
Field[3].ZMin = 0;
Field[3].ZMax = 1;

Field[4] = Min;
Field[4].FieldsList = {2, 3};

Background Field = 4;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

Mesh 3;