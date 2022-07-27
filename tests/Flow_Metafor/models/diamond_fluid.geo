/* Naca0012 */

// Geometry
DefineConstant[ xLgt = { 5.0, Name "Domain length (x-dir)" }  ];
DefineConstant[ yLgt = { 5.0, Name "Domain length (y-dir)" }  ];

// Mesh
DefineConstant[ msF = { 1.0, Name "Farfield mesh size" }  ];
DefineConstant[ msA = { 0.01, Name "Airfoil mesh size" }  ];

/**************
    Geometry  
 **************/

// Airfoil
Te = 1;
Le = 3;

Point(1) = {1,0,0,msA};
Point(2) = {.5,.04,0,msA};
Point(3) = {0,0,0,msA};
Point(4) = {.5,-0.04,0,msA};

Line(1) = {3,2};
Line(2) = {2,1};
Line(3) = {1,4};
Line(4) = {4,3};

// Farfield
Point(10001) = {1+xLgt, 0, 0,msF};
Point(10002) = {1+xLgt, yLgt, 0,msF};
Point(10003) = {-xLgt, yLgt, 0,msF};
Point(10004) = {-xLgt, 0, 0,msF};
Point(10005) = {-xLgt,-yLgt, 0,msF};
Point(10006) = {1+xLgt, -yLgt, 0,msF};

Line(10001) = {10001, 10002};
Line(10002) = {10002, 10003};
Line(10003) = {10003, 10004};
Line(10004) = {10004, 10005};
Line(10005) = {10005, 10006};
Line(10006) = {10006, 10001};

// Front and wake
Line(10007) = {Le, 10004};
Line(10008) = {Te, 10001};

// Internal field
Line Loop(10001) = {10008, 10001, 10002, 10003, -10007, 1, 2};
Line Loop(10002) = {10007, 10004, 10005, 10006, -10008, 3, 4};
Plane Surface(10001) = {10001};
Plane Surface(10002) = {10002};

/*************************
    Physical Groups
 *************************/

Physical Point("te") = {Te};
Physical Line("upstream") = {10003, 10004};
Physical Line("side") = {10002, 10005};
Physical Line("downstream") = {10001};
Physical Line("downstream_") = {10006};
Physical Line("airfoil") = {1, 2};
Physical Line("airfoil_") = {3, 4};
Physical Line("wake") = {10008};
Physical Surface("field") = {10001};
Physical Surface("field_") = {10002};
