// Gmsh project created on Fri Jan 09 15:36:17 2015
N = 40;
Ncurve = 6;
Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {0, 1, 0, 1.0};
Point(5) = {0.5, 0.195, 0, 1.0};
Point(6) = {0.505, 0.2, 0, 1.0};
Delete {
  Point{6, 5};
}
Point(5) = {0.5, 0.15, 0, 1.0};
Point(6) = {0.55, 0.2, 0, 1.0};
Point(7) = {0.5, 0.25, 0, 1.0};
Point(8) = {0.45, 0.2, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Circle(5) = {5, 6, 7};
Point(9) = {0.5, 0.2, 0, 1.0};
Delete {
  Line{5};
}
Circle(5) = {5, 9, 7};
Circle(6) = {7, 9, 5};
Delete {
  Point{8};
}
Delete {
  Point{6};
}
Physical Line("Wall",7) = {1, 2, 4};
Physical Line("Free",8) = {3};
Physical Line("Cylinder",9) = {6, 5};

Line Loop(10) = {1, 2, 3, 4};
Line Loop(11) = {5, 6};

Plane Surface(12) = {10, 11};
Plane Surface(13) = {11};

Physical Surface("Face1",14) = {12};
Physical Surface("Face2",15) = {13};

Transfinite Line {1} = N;
Transfinite Line {2} = N;
Transfinite Line {3} = N;
Transfinite Line {4} = N;
Transfinite Line {5} = Ncurve;
Transfinite Line {6} = Ncurve;
