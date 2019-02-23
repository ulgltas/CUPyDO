/* Diamond airfoil */

// Mesh
DefineConstant[ Nc = { 40, Name "Number of elements chordwise" }  ];

/**************
    Geometry  
 **************/

// Airfoil
Te = 1;
Le = 3;
Mc = 2;

Point(0) = {.5,0,0};
Point(1) = {1,0,0};
Point(2) = {.5,.04,0};
Point(3) = {0,0,0};
Point(4) = {.5,-0.04,0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line(5) = {Te,0};
Line(6) = {2,0};
Line(7) = {Le,0};
Line(8) = {4,0};

Line Loop(1) = {1,6,-5};
Surface(1) = {1};
Line Loop(2) = {2,7,-6};
Surface(2) = {2};
Line Loop(3) = {3,8,-7};
Surface(3) = {3};
Line Loop(4) = {4,5,-8};
Surface(4) = {4};

/***********
    Mesh
 ***********/
 
Transfinite Line{1,2,3,4,5,7} = Nc/2+1;
Transfinite Line{6,8} = 3;
Transfinite Surface{1,2,3,4};
Recombine Surface{1,2,3,4};

/*************************
    Physical Groups
 *************************/

Physical Point("te",101) = {Te};
Physical Point("le",103) = {Le};
Physical Point("mc",102) = {Mc};
Physical Line("airfoil",111) = {1,2,3,4};
Physical Surface("material",121) = {1,2,3,4};
