l = 0.2;  	//Tank half width
a = 0.02; 	//Cylinder Radius
h0 = 2*a; 	//Cylinder deepth 
h = h0+2*a+6*a; //Tank height without barrier
barrier = 2*a;  //Barrier height

N = 101;
Ncurve = (8*a/h)*N;

Point(1) = {-l, 0, 0, 1.0};
Point(2) = {l, 0, 0, 1.0};
Point(3) = {l, h, 0, 1.0};
Point(4) = {-l, h, 0, 1.0};
Point(5) = {0, h-(h0+2*a), 0, 1.0};
Point(6) = {0, h-h0, 0, 1.0};
Point(7) = {0, h-(h0+a), 0, 1.0};
Point(8) = {-l, h+barrier, 0, 1.0};
Point(9) = {l, h+barrier, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Circle(5) = {5, 7, 6};
Circle(6) = {6, 7, 5};
Line(7) = {4, 8};
Line(8) = {3, 9};

Line Loop(9) = {1, 2, 3, 4};
Line Loop(10) = {5, 6};
Plane Surface(11) = {9, 10};
Plane Surface(12) = {10};

Physical Line(13) = {1, 2, 4};
Physical Line(14) = {3};
Physical Line(15) = {7, 8};
Physical Line(16) = {6, 5};
Physical Surface(17) = {11};
Physical Surface(18) = {12};

Transfinite Line {1} = 1.5*(2*l/h)*N;
Transfinite Line {2} = 1.2*N;
Transfinite Line {3} = (2*l/h)*N;
Transfinite Line {4} = 1.2*N;
Transfinite Line {5} = Ncurve;
Transfinite Line {6} = Ncurve;
Transfinite Line {7} = 1.2*(barrier/h)*N+2;
Transfinite Line {8} = 1.2*(barrier/h)*N+2;