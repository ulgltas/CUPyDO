// Copyright 2018 University of Liège
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
//     http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License. 

// Gmsh project created on Wed Sep 14 18:09:29 2016

// --- Wing geometrical parameters ---
span = 0.762;
taper = 0.659;
AR = 1.65;
c_r = 0.559;		//root cord
c_t = taper*c_r;	//tip cord
sweep = Pi/4;
area = (c_t + c_r)*span/2.0;
delta = span*Tan(sweep) - c_t/4.0;	//parameter used to apply the sweep
DX = c_r/4.0+delta;			//another parameter used to apply the sweep

R_LE = 0.00102;		//Leading Edge radius
alpha = 115*Pi/180;	//opening angle of the leading edge

R_BL_r = 0.03;		//BL radius at root
R_BL_t = 0.02;		//BL radius at tip
alpha_2 = 115*Pi/180;   //opening angle of the leading edge BL
beta = 45*Pi/180;	//opening angle of the trailing edge BL

N_wallBL = 11;		//number of points in the BL (normal direction)
N_far = 11;
N_EI = 16;		//number of points on the extrados and intrados
N_LE = 2;		//number of point on the leading edge
N_TE = 2;		//number of points on the trailing edge
N_LEtip = 2;		//number of points on the LE radius
N_TEtip = N_LEtip;		//number of points on the TE radius
N_span = 32;		//number of points along the wing span

/*
q_LE = 0.6;
q_TE = 0.3;
q_span = 0.8;
*/


q_LE = 1.0;
q_TE = 1.0;
q_span = 1.0;




//Set of normally unused parameters, can be ignored...
finnesse_root = 0.05;
finnesse_tip = 0.002;
finnesse_LE = 0.002;
finnesse_TE = 0.002;
q_far = 0.7;

// --- Airfoil at the root (points cloud) ---
Point(1) = {0.000000, 0.000000, 0.000000, finnesse_LE};
Point(2) = {0.005000, 0.0, 0.003040, finnesse_LE};
Point(3) = {0.007500, 0.0, 0.003680, finnesse_LE};
Point(4) = {0.012500, 0.0, 0.004690, finnesse_LE};
Point(5) = {0.025000,0.0, 0.006470, finnesse_LE};
Point(6) = {0.050000, 0.0, 0.008750, finnesse_LE};
Point(7) = {0.075000, 0.0, 0.010590, finnesse_LE};
Point(8) = {0.100000, 0.0, 0.012130, finnesse_root};
Point(9) = {0.150000,0.0, 0.014590, finnesse_root};
Point(10) = {0.200000, 0.0, 0.016450, finnesse_root};
Point(11) = {0.250000, 0.0, 0.017880, finnesse_root};
Point(12) = {0.300000, 0.0, 0.018920, finnesse_root};
Point(13) = {0.350000, 0.0, 0.019620, finnesse_root};
Point(14) = {0.400000, 0.0, 0.019970, finnesse_root};
Point(15) = {0.450000, 0.0, 0.019960, finnesse_root};
Point(16) = {0.500000, 0.0, 0.019540, finnesse_root};
Point(17) = {0.550000, 0.0, 0.018680, finnesse_root};
Point(18) = {0.600000, 0.0,0.017430, finnesse_root};
Point(19) = {0.650000, 0.0, 0.015860, finnesse_root};
Point(20) = {0.700000, 0.0, 0.014020, finnesse_root};
Point(21) = {0.750000, 0.0, 0.011950, finnesse_root};
Point(22) = {0.800000, 0.0, 0.009670, finnesse_root};
Point(23) = {0.850000, 0.0, 0.007290, finnesse_root};
Point(24) = {0.900000, 0.0, 0.004900, finnesse_root};
Point(25) = {0.950000, 0.0, 0.002500, finnesse_root};
Point(26) = {1.000000, 0.0, 0.000090, finnesse_root};
Point(27) = {1.000000, 0.0, -0.000090, finnesse_root};
//Point(27) = {1.000000, 0.0, 0.000000, finnesse_root};
Point(28) = {0.950000, 0.0, -0.002500, finnesse_root};
Point(29) = {0.900000, 0.0, -0.004900, finnesse_root};
Point(30) = {0.850000, 0.0, -0.007290, finnesse_root};
Point(31) = {0.800000, 0.0, -0.009670, finnesse_root};
Point(32) = {0.750000, 0.0, -0.011950, finnesse_root};
Point(33) = {0.700000, 0.0, -0.014020, finnesse_root};
Point(34) = {0.650000, 0.0, -0.015860, finnesse_root};
Point(35) = {0.600000, 0.0, -0.017430, finnesse_root};
Point(36) = {0.550000,0.0, -0.018680, finnesse_root};
Point(37) = {0.500000, 0.0, -0.019540, finnesse_root};
Point(38) = {0.450000, 0.0, -0.019960, finnesse_root};
Point(39) = {0.400000, 0.0, -0.019970, finnesse_root};
Point(40) = {0.350000, 0.0, -0.019620, finnesse_root};
Point(41) = {0.300000, 0.0, -0.018920, finnesse_root};
Point(42) = {0.250000, 0.0, -0.017880, finnesse_root};
Point(43) = {0.200000, 0.0, -0.016450, finnesse_root};
Point(44) = {0.150000, 0.0, -0.014590, finnesse_root};
Point(45) = {0.100000, 0.0, -0.012130, finnesse_root};
Point(46) = {0.075000, 0.0, -0.010590, finnesse_LE};
Point(47) = {0.050000, 0.0, -0.008750, finnesse_LE};
Point(48) = {0.025000, 0.0, -0.006470, finnesse_LE};
Point(49) = {0.012500, 0.0, -0.004690, finnesse_LE};
Point(50) = {0.007500, 0.0, -0.003680, finnesse_LE};
Point(51) = {0.005000, 0.0, -0.003040, finnesse_LE};
Point(52) = {0.00102, 0.0, 0.000000, finnesse_root};
Point(53) = {R_LE*Cos(alpha)+R_LE, 0.0, R_LE*Sin(alpha), finnesse_LE};
Point(54) = {R_LE*Cos(alpha)+R_LE, 0.0, -R_LE*Sin(alpha), finnesse_LE};
Point(109) = {1+0.000090, 0.0, 0.0, finnesse_TE};
Point(110) = {1, 0.0, 0.0, finnesse_TE};

// --- Airfoil at the tip (points cloud) ---
Point(55) = {0.000000, span, 0.000000, finnesse_tip};
Point(56) = {0.005000, span, 0.003040, finnesse_tip};
Point(57) = {0.007500, span, 0.003680, finnesse_tip};
Point(58) = {0.012500, span, 0.004690, finnesse_tip};
Point(59) = {0.025000,span, 0.006470, finnesse_tip};
Point(60) = {0.050000, span, 0.008750, finnesse_tip};
Point(61) = {0.075000, span, 0.010590, finnesse_tip};
Point(62) = {0.100000, span, 0.012130, finnesse_tip};
Point(63) = {0.150000,span, 0.014590, finnesse_tip};
Point(64) = {0.200000, span, 0.016450, finnesse_tip};
Point(65) = {0.250000, span, 0.017880, finnesse_tip};
Point(66) = {0.300000, span, 0.018920, finnesse_tip};
Point(67) = {0.350000, span, 0.019620, finnesse_tip};
Point(68) = {0.400000, span, 0.019970, finnesse_tip};
Point(69) = {0.450000, span, 0.019960, finnesse_tip};
Point(70) = {0.500000, span, 0.019540, finnesse_tip};
Point(71) = {0.550000, span, 0.018680, finnesse_tip};
Point(72) = {0.600000, span,0.017430, finnesse_tip};
Point(73) = {0.650000, span, 0.015860, finnesse_tip};
Point(74) = {0.700000, span, 0.014020, finnesse_tip};
Point(75) = {0.750000, span, 0.011950, finnesse_tip};
Point(76) = {0.800000, span, 0.009670, finnesse_tip};
Point(77) = {0.850000, span, 0.007290, finnesse_tip};
Point(78) = {0.900000, span, 0.004900, finnesse_tip};
Point(79) = {0.950000, span, 0.002500, finnesse_tip};
Point(80) = {1.000000, span, 0.000090, finnesse_tip};
Point(81) = {1.000000, span, -0.000090, finnesse_tip};
//Point(81) = {1.000000, span, 0.000000, finnesse_tip};
Point(82) = {0.950000, span, -0.002500, finnesse_tip};
Point(83) = {0.900000, span, -0.004900, finnesse_tip};
Point(84) = {0.850000, span, -0.007290, finnesse_tip};
Point(85) = {0.800000, span, -0.009670, finnesse_tip};
Point(86) = {0.750000, span, -0.011950, finnesse_tip};
Point(87) = {0.700000, span, -0.014020, finnesse_tip};
Point(88) = {0.650000, span, -0.015860, finnesse_tip};
Point(89) = {0.600000, span, -0.017430, finnesse_tip};
Point(90) = {0.550000,span, -0.018680, finnesse_tip};
Point(91) = {0.500000, span, -0.019540, finnesse_tip};
Point(92) = {0.450000, span, -0.019960, finnesse_tip};
Point(93) = {0.400000, span, -0.019970, finnesse_tip};
Point(94) = {0.350000, span, -0.019620, finnesse_tip};
Point(95) = {0.300000, span, -0.018920, finnesse_tip};
Point(96) = {0.250000, span, -0.017880, finnesse_tip};
Point(97) = {0.200000, span, -0.016450, finnesse_tip};
Point(98) = {0.150000, span, -0.014590, finnesse_tip};
Point(99) = {0.100000, span, -0.012130, finnesse_tip};
Point(100) = {0.075000, span, -0.010590, finnesse_tip};
Point(101) = {0.050000, span, -0.008750, finnesse_tip};
Point(102) = {0.025000, span, -0.006470, finnesse_tip};
Point(103) = {0.012500, span, -0.004690, finnesse_tip};
Point(104) = {0.007500, span, -0.003680, finnesse_tip};
Point(105) = {0.005000, span, -0.003040, finnesse_tip};
Point(106) = {0.00102, span, 0.000000, finnesse_tip};
Point(107) = {R_LE*Cos(alpha)+R_LE, span, R_LE*Sin(alpha), finnesse_tip};
Point(108) = {R_LE*Cos(alpha)+R_LE, span, -R_LE*Sin(alpha), finnesse_tip};
Point(111) = {1+0.000090, span, 0.0, finnesse_TE};
Point(112) = {1, span, 0.0, finnesse_TE};

// --- Airfoil at the root with splines ---
Spline(1) = {7, 6, 5, 4, 3, 2, 53};
Circle(2) = {53, 52, 1};
Circle(3) = {1, 52, 54};
Spline(4) = {54, 51, 50, 49, 48, 47, 46};
Spline(5) = {7:25};
Spline(6) = {28:46};
Line(7) = {25, 26};
Line(8) = {28, 27};
Circle(9) = {26, 110, 109};
Circle(10) = {109, 110, 27};


// --- Airfoil at the tip with splines ---
Spline(11) = {107, 56, 57, 58, 59, 60, 61};
Spline(12) = {108, 105, 104, 103, 102, 101, 100};
Circle(13) = {107, 106, 55};
Circle(14) = {55, 106, 108};
Spline(15) = {61:79};
Spline(16) = {82:100};
Line(17) = {79, 80};
Line(18) = {82, 81};
Circle(19) = {80, 112, 111};
Circle(20) = {111, 112, 81};

// --- Rescale root and tip : apply taper ratio ---
Dilate {{0, 0, 0}, c_r} {
  Line{5, 1, 2, 3, 4, 6, 8, 7, 9, 10};
}
Dilate {{0, span, 0}, c_t} {
  Line{14, 17, 19, 11, 12, 13, 16, 15, 18, 20};
}

// --- Mesh the root airfoil ---
//Geometry
Point(113) = {0.041925, 0.0, 0.0, 1.0};
Point(114) = {0.53105, 0.0, 0.0, 1.0};
Line(21) = {7, 113};
Line(22) = {113, 1};
Line Loop(23) = {1, 2, -22, -21};
Plane Surface(24) = {23};
Line(25) = {46, 113};
Line Loop(26) = {22, 3, 4, 25};
Plane Surface(27) = {26};
Line(28) = {25, 114};
Line(29) = {28, 114};
Line(30) = {114, 109};
Line(31) = {114, 113};
Line Loop(32) = {5, 28, 31, -21};
Plane Surface(33) = {32};
Line Loop(34) = {31, -25, -6, 29};
Plane Surface(35) = {34};
Line Loop(36) = {7, 9, -30, -28};
Plane Surface(37) = {36};
Line Loop(38) = {30, 10, -8, 29};
Plane Surface(39) = {38};
//Discretisation
Transfinite Line {21, 2, 25, 3} = N_LEtip Using Progression 1;
Transfinite Line {1, 22} = N_LE Using Progression q_LE;
Transfinite Line {4} = N_LE Using Progression 1.0/q_LE;
Transfinite Line {9, 10, 28, 29} = N_TEtip Using Progression 1;
Transfinite Line {8, 7, 30} = N_TE Using Progression q_TE;
Transfinite Line {5, 31, 6} = N_EI Using Progression 1;
Transfinite Surface {24} = {7, 113, 1, 53};
Transfinite Surface {27} = {113, 46, 54, 1};
Transfinite Surface {33} = {113, 7, 25, 114};
Transfinite Surface {35} = {114, 28, 46, 113};
Transfinite Surface {37} = {114, 25, 26, 109};
Transfinite Surface {39} = {109, 27, 28, 114};


// --- Mesh the tip airfoil ---
//Geometry
Point(115) = {0.027628575, span, 0.0, 1.0};
Point(116) = {0.34996195, span, 0.0, 1.0};
Line(50) = {100, 115};
Line(51) = {61, 115};
Line(52) = {115, 55};
Line(53) = {79, 116};
Line(54) = {82, 116};
Line(55) = {116, 115};
Line(56) = {116, 111};
//Apply the sweep back
Translate {DX, 0, 0} {
  Line{15, 11, 13, 14, 12, 52, 50, 51, 16, 55, 17, 18, 54, 53, 19, 20, 56};
}
Point(117) = {0.809871700391, 0.7618, 0, 1.0};
Point(118) = {0.809871700391, 0.76217, 0.0, 1.0};
Point(119) = {0.809871700391, 0.762, 0.0, 1.0};
Ellipse(57) = {118, 119, 107, 107};
Ellipse(58) = {55, 117, 55, 118};
Ellipse(59) = {108, 119, 108, 118};
Extrude {{1, 0, 0}, {0, 0.762, 0}, -Pi/2} {
  Line{15};
}
Extrude {{1, 0, 0}, {0, 0.762, 0}, Pi/2} {
  Line{16};
}
Rotate {{1, 0, 0}, {0, 0.762, 0}, -Pi/2} {
  Duplicata { Point{60, 59, 58, 57, 56}; }
}
Spline(68) = {118, 143, 142, 141, 140, 139, 120};
Line Loop(69) = {61, -68, 57, 11};
Ruled Surface(70) = {69};
Line Loop(71) = {66, -68, -59, 12};
Ruled Surface(72) = {71};
Line Loop(73) = {57, 13, 58};
Ruled Surface(74) = {73};
Line Loop(75) = {59, -58, 14};
Ruled Surface(76) = {75};
Compound Surface(77) = {74, 76};
Extrude {{1, 0, 0}, {0, 0.762, 0}, -Pi/2} {
  Line{17, 19};
}
Extrude {{1, 0, 0}, {0, 0.762, 0}, Pi/2} {
  Line{18, 20};
}
Compound Surface(92) = {84, 91};
Line Loop(93) = {53, -54, 65, -62};
Plane Surface(94) = {93};
Line Loop(95) = {53, 56, -19, -17};
Plane Surface(96) = {95};
Line Loop(97) = {54, 56, 20, -18};
Plane Surface(98) = {97};
Surface Loop(99) = {94, 96, 98, 88, 81, 92};
Volume(100) = {99};
Line Loop(101) = {51, -50, 66, -61};
Plane Surface(102) = {101};
Line Loop(103) = {51, -55, -53, -15};
Plane Surface(104) = {103};
Line Loop(105) = {55, -50, -16, 54};
Plane Surface(106) = {105};
Surface Loop(107) = {102, 104, 106, 67, 63, 94};
Volume(108) = {107};
Line Loop(109) = {11, 51, 52, -13};
Plane Surface(110) = {109};
Line Loop(111) = {52, 14, 12, 50};
Plane Surface(112) = {111};
Surface Loop(113) = {112, 110, 70, 72, 102, 77};
Volume(114) = {113};

//Discretisation
Transfinite Line {19, 20, 80, 87} = N_TEtip Using Progression 1;
Transfinite Line {17, 78, 18, 56} = N_TE Using Progression q_TE;
Transfinite Line {53, 62, 65, 54} = N_TEtip Using Progression 1;
Transfinite Line {16, 60, 15, 55} = N_EI Using Progression 1;
Transfinite Line {61, 66, 50, 51} = N_TEtip Using Progression 1;
Transfinite Line {52} = N_LE Using Progression q_LE;
Transfinite Line {11, 68, 12} = N_LE Using Progression 1.0/q_LE;
Transfinite Line {57, 13, 14, 59} = N_TEtip Using Progression 1;
Transfinite Surface {92} = {80, 145, 81, 111};
Transfinite Surface {81} = {79, 138, 145, 80};
Transfinite Surface {88} = {81, 145, 138, 82};
Transfinite Surface {96} = {116, 79, 80, 111};
Transfinite Surface {98} = {111, 81, 82, 116};
Transfinite Surface {94} = {79, 116, 82, 138};
Transfinite Surface {63} = {61, 120, 138, 79};
Transfinite Surface {67} = {138, 82, 100, 120};
Transfinite Surface {102} = {61, 120, 100, 115};
Transfinite Surface {104} = {115, 61, 79, 116};
Transfinite Surface {106} = {116, 82, 100, 115};
Transfinite Surface {70} = {120, 61, 107, 118};
Transfinite Surface {72} = {118, 108, 100, 120};
Transfinite Surface {77} = {107, 55, 108, 118};
Transfinite Surface {110} = {61, 115, 55, 107};
Transfinite Surface {112} = {55, 108, 100, 115};
Transfinite Volume{100} = {79, 116, 82, 138, 80, 111, 81, 145};
Transfinite Volume{108} = {79, 116, 82, 138, 61, 115, 100, 120};
Transfinite Volume{114} = {107, 55, 108, 118, 61, 115, 100, 120};


// --- Wing skin and internal mesh ---
//Geometry
Line(118) = {7, 61};
Line(119) = {25, 79};
Line(120) = {109, 111};
Line(121) = {114, 116};
Line(122) = {113, 115};
Line(123) = {46, 100};
Line(124) = {28, 82};
Line(125) = {1, 55};
Line(126) = {26, 80};
Line(127) = {27, 81};
Line Loop(128) = {10, 127, -20, -120};
Ruled Surface(129) = {128};
Line Loop(130) = {127, -18, -124, 8};
Plane Surface(131) = {130};

Line Loop(132) = {121, 56, -120, -30};
Plane Surface(133) = {132};
Line Loop(134) = {124, 54, -121, -29};
Plane Surface(135) = {134};
Surface Loop(136) = {135, 131, 129, 39, 133, 98};
Volume(137) = {136};
Line Loop(138) = {9, 120, -19, -126};
Ruled Surface(139) = {138};
Line Loop(140) = {7, 126, -17, -119};
Plane Surface(141) = {140};
Line Loop(142) = {119, 53, -121, -28};
Plane Surface(143) = {142};
Surface Loop(144) = {37, 141, 139, 143, 133, 96};
Volume(145) = {144};

Line Loop(146) = {123, -16, -124, 6};
Ruled Surface(147) = {146};
Line Loop(148) = {31, 122, -55, -121};
Plane Surface(149) = {148};
Line Loop(150) = {123, 50, -122, -25};
Plane Surface(151) = {150};
Surface Loop(152) = {35, 147, 151, 149, 135, 106};
Volume(153) = {152};

Line Loop(154) = {119, -15, -118, 5};
Ruled Surface(155) = {154};
Line Loop(156) = {118, 51, -122, -21};
Plane Surface(157) = {156};
Surface Loop(158) = {33, 155, 157, 143, 149, 104};
Volume(159) = {158};

Line(160) = {53, 107};
Line(161) = {54, 108};
Line Loop(162) = {3, 161, -14, -125};
Ruled Surface(163) = {162};
Line Loop(164) = {4, 123, -12, -161};
Ruled Surface(165) = {164};
Line Loop(166) = {122, 52, -125, -22};
Plane Surface(167) = {166};
Surface Loop(168) = {27, 163, 165, 167, 151, 112};
Volume(169) = {168};
Line Loop(170) = {160, 13, -125, -2};
Ruled Surface(171) = {170};
Line Loop(172) = {1, 160, 11, -118};
Ruled Surface(173) = {172};
Surface Loop(174) = {24, 173, 171, 157, 167, 110};
Volume(175) = {174};

//Discretisation
Transfinite Line {125, 123, 118, 122, 124, 121, 119, 120, 126, 127, 160, 161} = N_span Using Progression q_span;
Transfinite Surface {129} = {109, 27, 81, 111};
Transfinite Surface {131} = {27, 28, 82, 81};
Transfinite Surface {133} = {109, 114, 116, 111};
Transfinite Surface {135} = {114, 28, 82, 116};
Transfinite Surface {139} = {26, 109, 111, 80};
Transfinite Surface {141} = {26, 25, 79, 80};
Transfinite Surface {143} = {25, 114, 116, 79};
Transfinite Surface {147} = {28, 46, 100, 82};
Transfinite Surface {149} = {114, 113, 115, 116};
Transfinite Surface {151} = {100, 115, 113, 46};
Transfinite Surface {155} = {25, 7, 61, 79};
Transfinite Surface {157} = {61, 115, 113, 7};
Transfinite Surface {163} = {1, 54, 108, 55};
Transfinite Surface {165} = {46, 54, 108, 100};
Transfinite Surface {167} = {113, 1, 55, 115};
Transfinite Surface {171} = {53, 1, 55, 107};
Transfinite Surface {173} = {7, 53, 107, 61};
Transfinite Volume{137} = {114, 28, 27, 109, 116, 82, 81, 111};
Transfinite Volume{145} = {25, 26, 109, 114, 79, 80, 111, 116};
Transfinite Volume{153} = {113, 114, 28, 46, 115, 116, 82, 100};
Transfinite Volume{159} = {7, 113, 114, 25, 61, 115, 116, 79};
Transfinite Volume{169} = {113, 46, 54, 1, 115, 100, 108, 55};
Transfinite Volume{175} = {7, 113, 1, 53, 61, 115, 55, 107};

Mesh.RecombineAll = 1;

/*
//Physical groups
Physical Surface(176) = {24, 27, 33, 35, 37, 39};	//clamped face
Physical Surface(177) = {173, 171, 163, 155, 141, 139, 129, 131, 165, 147, 71, 63, 70, 91, 95, 103, 110, 114};						//free surface

Physical Line(178) = {15, 17, 11};			//traction line
Physical Volume(179) = {145, 137, 159, 153, 175, 169, 79, 99, 116};	//internal mesh
Physical Point(180) = {70};				//extract the displacement
*/


Physical Surface(176) = {37, 39, 33, 35, 24, 27};  	//clamped face
Physical Surface(177) = {173, 155, 141, 171, 163, 77, 70, 72, 63, 67, 81, 88, 92, 165, 147, 131, 139, 129};							//free surface
Physical Line(178) = {60};				//traction line
Physical Volume(179) = {169, 175, 159, 153, 137, 145, 114, 108, 100};	//internal mesh
Physical Point(180) = {55};				//extract the displacement at LE tip
Physical Point(181) = {111};				//extract the displacement at TE tip
Physical Surface(182) = {147, 165, 131};		//traction surface
Physical Line(183) = {125};                 //Leading edge
Physical Line(184) = {120};                 //Trailing edge
