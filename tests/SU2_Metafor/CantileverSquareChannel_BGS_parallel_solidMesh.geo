// Gmsh project created on Wed Jan 13 14:20:12 2016
//Parameter
H = 0.01;
L = 0.04;
t = 0.0006;
U = 0.513;
rho = 1.18e-4;
mu = 1.82e-5;

//Body definition
Point(5) = {H/2, t/2, 0.0,1.0};
Point(6) = {H/2, -t/2, 0.0, 1.0};
Point(7) = {H/2+L, t/2, 0.0,1.0};
Point(8) = {H/2+L, -t/2, 0.0, 1.0};

//Body boundaries
Line(5) = {6, 8};
Line(6) = {8, 7};
Line(7) = {7, 5};
Line(65) = {5, 6};
Line Loop(66) = {5, 6, 7, 65};
Plane Surface(67) = {66};

//Internal domain discretisation
Transfinite Line {7, 5} = 241 Using Progression 1;
Transfinite Line {6, 65} = 11 Using Progression 1;
Transfinite Surface {67};
Recombine Surface {67}; // tri => quads

//Physical boundaries and physical domain
Physical Line(101) = {65};         // clamped side of the beam
Physical Line(102) = {7, 5, 6};    // free surface of the beam
Physical Line(103) = {7};          // upper surface of the beam (for tests only)
Physical Surface(100) = {67};      // meshed beam
Physical Point(104) = {7};         // upper right vertex
