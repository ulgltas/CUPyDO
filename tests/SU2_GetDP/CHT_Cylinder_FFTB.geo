// Gmsh project created on Wed Feb 22 13:55:24 2017

D = 1.0/1;
Dc = 0.5*D;

N_circ = 31;
N_rad = 21;

//Center
Point(1) = {0, 0, 0, 1.0};

//Inner diameter
Point(2) = {Dc/2, 0, 0, 1.0};
Point(3) = {0, Dc/2, 0, 1.0};
Point(4) = {-Dc/2, 0, 0, 1.0};
Point(5) = {0, -Dc/2, 0, 1.0};

//Outer diamter
Point(6) = {D/2, 0, 0, 1.0};
Point(7) = {0, D/2, 0, 1.0};
Point(8) = {-D/2, 0, 0, 1.0};
Point(9) = {0, -D/2, 0, 1.0};

//Domain definition
Circle(1) = {3, 1, 2};
Circle(2) = {2, 1, 5};
Circle(3) = {5, 1, 4};
Circle(4) = {4, 1, 3};
Circle(5) = {7, 1, 6};
Circle(6) = {6, 1, 9};
Circle(7) = {9, 1, 8};
Circle(8) = {8, 1, 7};
Line(9) = {7, 3};
Line(10) = {6, 2};
Line(11) = {9, 5};
Line(12) = {8, 4};
Line Loop(13) = {8, 9, -4, -12};
Plane Surface(14) = {13};
Line Loop(15) = {9, 1, -10, -5};
Plane Surface(16) = {15};
Line Loop(17) = {10, 2, -11, -6};
Plane Surface(18) = {17};
Line Loop(19) = {3, -12, -7, 11};
Plane Surface(20) = {19};

//Domain discretization
Transfinite Line {8, 4, 7, 6, 5, 1, 2, 3} = N_circ Using Progression 1;
Transfinite Line {12, 9, 10, 11} = N_rad Using Progression 1;
Transfinite Surface {14} = {8, 4, 3, 7};
Transfinite Surface {16} = {7, 3, 2, 6};
Transfinite Surface {18} = {2, 6, 9, 5};
Transfinite Surface {20} = {4, 8, 9, 5};
Recombine Surface {14, 20, 18, 16};

//Physical groups
Physical Line("inner", 101) = {4, 3, 2, 1};  //inner
Physical Line("outer", 102) = {8, 7, 6, 5}; //outer
Physical Surface("internal", 103) = {14, 16, 18, 20};   //internal
Physical Point(104) = {3};      //inner point
Physical Point(105) = {7};      //outer point
Physical Line(106) = {9};       //extractor through the thickness
