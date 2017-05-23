// Gmsh project created on Wed Jan 13 14:20:12 2016
//Parameter
H = 0.01;
L = 0.04;
t = 0.0006;
U = 0.513;
rho = 1.18e-4;
mu = 1.82e-5;

//Body definition
Point(1) = {H/2, H/2, 0.0,1.0};
Point(2) = {-H/2, H/2, 0.0,1.0};
Point(3) = {-H/2, -H/2, 0.0,1.0};
Point(4) = {H/2, -H/2, 0.0,1.0};
Point(5) = {H/2, t/2, 0.0,1.0};
Point(6) = {H/2, -t/2, 0.0, 1.0};
Point(7) = {H/2+L, t/2, 0.0,1.0};
Point(8) = {H/2+L, -t/2, 0.0, 1.0};

//External domain
Point(9) = {14.5*H, 6*H, 0.0,1.0};
Point(10) = {14.5*H, -6*H, 0.0, 1.0};
Point(11) = {-5*H, 6*H, 0.0,1.0};
Point(12) = {-5*H, -6*H, 0.0, 1.0};
Point(13) = {5*H, 6*H, 0.0,1.0};
Point(14) = {5*H, -6*H, 0.0, 1.0};

//Body boundaries
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 6};
Line(5) = {6, 8};
Line(6) = {8, 7};
Line(7) = {7, 5};
Line(8) = {5, 1};

//External boundaries
Line(9) = {11, 12};
Line(10) = {12, 14};
Line(11) = {14, 10};
Line(12) = {10, 9};
Line(13) = {9, 13};
Line(14) = {13, 11};

//Body extension definition
growth = 0.002;
Point(15) = {H/2+growth, H/2+growth, 0.0,1.0};
Point(16) = {-H/2-growth, H/2+growth, 0.0,1.0};
Point(17) = {-H/2-growth, -H/2-growth, 0.0,1.0};
Point(18) = {H/2+growth, -H/2-growth, 0.0,1.0};
Point(19) = {H/2 +growth, t/2+growth, 0.0,1.0};
Point(20) = {H/2+growth, -t/2-growth, 0.0, 1.0};
Point(21) = {(H/2+L)+growth, t/2+growth, 0.0,1.0};
Point(22) = {(H/2+L)+growth, -t/2-growth, 0.0, 1.0};

//Body extension boundaries
Line(15) = {15, 16};
Line(16) = {16, 17};
Line(17) = {17, 18};
Line(18) = {18, 20};
Line(19) = {20, 22};
Line(20) = {22, 21};
Line(21) = {21, 19};
Line(22) = {19, 15};

//Internal connexion lines
Line(23) = {16, 2};
Line(24) = {17, 3};
Line(25) = {18, 4};
Line(26) = {15, 1};
Line(27) = {19, 5};
Line(28) = {20, 6};
Line(29) = {21, 7};
Line(30) = {22, 8};

//external connexion lines
Line(31) = {11, 16};
Line(32) = {13, 15};
Line(33) = {9, 21};
Line(34) = {10, 22};
Line(35) = {14, 18};
Line(36) = {12, 17};

//External domain definition
Line Loop(37) = {9, 36, -16, -31};
Plane Surface(38) = {37};
Line Loop(39) = {31, -15, -32, 14};
Plane Surface(40) = {39};
Line Loop(41) = {13, 32, -22, -21, -33};
Plane Surface(42) = {41};
Line Loop(43) = {12, 33, -20, -34};
Plane Surface(44) = {43};
Line Loop(45) = {11, 34, -19, -18, -35};
Plane Surface(46) = {45};
Line Loop(47) = {35, -17, -36, 10};
Plane Surface(48) = {47};

//Internal domain definition
Line Loop(49) = {16, 24, -2, -23};
Plane Surface(50) = {49};
Line Loop(51) = {15, 23, -1, -26};
Plane Surface(52) = {51};
Line Loop(53) = {26, -8, -27, 22};
Plane Surface(54) = {53};
Line Loop(55) = {21, 27, -7, -29};
Plane Surface(56) = {55};
Line Loop(57) = {29, -6, -30, 20};
Plane Surface(58) = {57};
Line Loop(59) = {19, 30, -5, -28};
Plane Surface(60) = {59};
Line Loop(61) = {4, -28, -18, 25};
Plane Surface(62) = {61};
Line Loop(63) = {25, -3, -24, 17};
Plane Surface(64) = {63};

//Internal domain discretisation
Transfinite Line {1, 2, 3, 15, 16, 17} = 31 Using Progression 1;
Transfinite Line {8, 4, 18, 22} = 16 Using Progression 1;
Transfinite Line {7, 5, 21, 19} = 121 Using Progression 1;
Transfinite Line {6, 20} = 6 Using Progression 1;
Transfinite Line {23, 24, 25, 28, 27, 26, 29, 30} = 11 Using Progression 1;
Transfinite Surface {50} = {17, 3, 2, 16};
Transfinite Surface {52} = {2, 1, 15, 16};
Transfinite Surface {54} = {15, 1, 5, 19};
Transfinite Surface {56} = {19, 5, 7, 21};
Transfinite Surface {58} = {21, 7, 8, 22};
Transfinite Surface {60} = {22, 8, 6, 20};
Transfinite Surface {62} = {4, 18, 20, 6};
Transfinite Surface {64} = {17, 18, 4, 3};
Recombine Surface {50, 52, 54, 56, 58, 60, 62, 64};

//External domain discretisation
Transfinite Line {9, 14, 13, 12, 11, 10} = 31 Using Progression 1;
Transfinite Line {31, 36, 35, 34, 33, 32} = 31 Using Progression 0.9;
Recombine Surface {38, 40, 42, 44, 46, 48};

//Physical boundaries and physical domain
/*Physical Line("inlet") = {9};
Physical Line("outlet") = {12};
Physical Line("slipWall") = {13, 14, 10, 11};
Physical Line("bodyWall") = {1, 2, 3, 4, 8, 7, 5, 6};
Physical Surface("internal") = {38, 40, 48, 42, 44, 46, 50, 52, 54, 64, 62, 56, 60, 58};*/ 
// [RB] metafor is not able to read strings as id yet
// [RB] I'll do that later
Physical Line(1) = {9};
Physical Line(2) = {12};
Physical Line(3) = {13, 14, 10, 11};
Physical Line(4) = {1, 2, 3, 4, 8, 7, 5, 6};
Physical Surface(5) = {38, 40, 48, 42, 44, 46, 50, 52, 54, 64, 62, 56, 60, 58};


// --- added by RB

// solid domain
Line(65) = {5, 6};
Transfinite Line {65} = 6 Using Progression 1;
Line Loop(66) = {5, 6, 7, 65};
Plane Surface(67) = {66};
Transfinite Surface {67};
Recombine Surface {67}; // tri => quads

Physical Line(101) = {65};         // clamped side of the beam
Physical Line(102) = {7, 5, 6};    // free surface of the beam
Physical Line(103) = {7};          // upper surface of the beam (for tests only)
Physical Surface(100) = {67};      // meshed beam

Physical Point(104) = {7};         // upper right vertex


