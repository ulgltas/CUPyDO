// Gmsh project created on Thu Apr 19 14:52:53 2018
N = 100;
Nside = N/10+1;

Point(1) = {-0.5, 0, 0, 1.0};
Point(2) = {0.5, 0, 0, 1.0};
Point(3) = {0.5, 1, 0, 1.0};
Point(4) = {-0.5, 1, 0, 1.0};
Point(5) = {-0.05, 0.1, 0, 1.0};
Point(6) = {0.05, 0.1, 0, 1.0};
Point(7) = {0.05, 0.2, 0, 1.0};
Point(8) = {-0.05, 0.2, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

Line Loop(1) = {1, 2, 3, 4};
Line Loop(2) = {8, 5, 6, 7};

Plane Surface(1) = {1, 2};

Physical Line(1) = {1, 2, 4};
Physical Line(2) = {3};
Physical Line(3) = {5, 6, 7, 8};
Physical Surface(4) = {1};

Transfinite Line {1} = N;
Transfinite Line {2} = N;
Transfinite Line {3} = N;
Transfinite Line {4} = N;
Transfinite Line {5} = Nside;
Transfinite Line {6} = Nside;
Transfinite Line {7} = Nside;
Transfinite Line {8} = Nside;
