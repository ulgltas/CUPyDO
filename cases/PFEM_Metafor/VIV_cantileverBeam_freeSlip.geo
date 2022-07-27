// Gmsh project created on Tue Mar 28 10:37:47 2017
//Parameter
H = 0.01;
L = 0.04;
t = 0.0006;
U = 0.513;
rho = 1.18e-4;
mu = 1.82e-5;
N = 11;

Point(1) = {-5*H, -6*H, 0, 1.0};
Point(2) = {14.5*H, -6*H, 0, 1.0};
Point(3) = {14.5*H, 6*H, 0, 1.0};
Point(4) = {-5*H, 6*H, 0, 1.0};
Point(5) = {0.5*H, -0.5*H, 0, 1.0};
Point(6) = {-0.5*H, -0.5*H, 0, 1.0};
Point(7) = {-0.5*H, 0.5*H, 0, 1.0};
Point(8) = {0.5*H, 0.5*H, 0, 1.0};
Point(9) = {0.5*H, 0.5*t, 0, 1.0};
Point(10) = {0.5*H, -0.5*t, 0, 1.0};
Point(11) = {0.5*H+L, -0.5*t, 0, 1.0};
Point(12) = {0.5*H+L, 0.5*t, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {10, 5};
Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 8};
Line(9) = {8, 9};
Line(10) = {9, 10};
Line(11) = {10, 11};
Line(12) = {11, 12};
Line(13) = {12, 9};

Line Loop(14) = {4, 1, 2, 3};
Line Loop(15) = {5, 6, 7, 8, 9, 10};
Line Loop(16) = {11, 12, 13, 10};
Plane Surface(17) = {14,15,-16};
Plane Surface(18) = {16};

Physical Line(19) = {1, 3};           //Upper and lower external boundaries
Physical Line(20) = {2};              //Outlet
Physical Line(21) = {4};     	      //Inlet
Physical Line(22) = {5, 6, 7, 8, 9};  //Square surface
Physical Line(23) = {11, 12, 13};     //Beam surface: FSI interface
Physical Line(24) = {10};             //Beam clamping
Physical Surface(25) = {17};   	      //Air
Physical Surface(26) = {18};          //Beam
Physical Surface(100) = {};           //Solid-solid contact

Transfinite Line{1} = 20*N;
Transfinite Line{2} = 20*N;
Transfinite Line{3} = 20*N;
Transfinite Line{4} = 10*N;
Transfinite Line{5} = ((H/2 - t/2)/H)*N+1;
Transfinite Line{6} = N;
Transfinite Line{7} = N;
Transfinite Line{8} = N;
Transfinite Line{9} = ((H/2 - t/2)/H)*N+1;
Transfinite Line{10} = 3;
Transfinite Line{11} = 5*N;
Transfinite Line{12} = 3;
Transfinite Line{13} = 5*N;

Transfinite Surface{18} = {9,10,11,12};
Recombine Surface{18};Rotate {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Surface{18, 17};
}
Rotate {{0, 0, 1}, {0, 0, 0}, -Pi/2} {
  Surface{17, 18};
}
Rotate {{0, 0, 1}, {0, 0, 0}, -Pi/2} {
  Surface{18, 17};
}
