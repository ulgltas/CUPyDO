/* AGARD 445 wing */
// Initially generated by unsRgridWingGen.m
// For gmsh 4, use line 220 instead of line 221 (Bezier <- BSpline)

// Parameters
// domain and mesh
DefineConstant[ xO = { -3, Name "Box origin (x)" }  ];
DefineConstant[ zO = { -3, Name "Box origin (z)" }  ];
DefineConstant[ xL = { 7, Name "Box length" }  ];
DefineConstant[ yL = { 3, Name "Box width" }  ];
DefineConstant[ zL = { 6, Name "Box height" }  ];
DefineConstant[ msLeRt = { 0.005, Name "Root leading edge mesh size" }  ];
DefineConstant[ msTeRt = { 0.005, Name "Root trailing edge mesh size" }  ];
DefineConstant[ msLeTp = { 0.003, Name "Tip leading edge mesh size" }  ];
DefineConstant[ msTeTp = { 0.003, Name "Tip trailing edge mesh size" }  ];
DefineConstant[ msF = { 1.0, Name "Farfield mesh size" }  ];

//// GEOMETRY


/// Points

// Airfoil 1: agard445, 51 points 
Point(1) = {0.559000,0.000000,0.000000,msTeRt};
Point(2) = {0.531050,0.000000,0.001398};
Point(3) = {0.503100,0.000000,0.002739,msTeRt};
Point(4) = {0.475150,0.000000,0.004075};
Point(5) = {0.447200,0.000000,0.005406};
Point(6) = {0.419250,0.000000,0.006680};
Point(7) = {0.391300,0.000000,0.007837};
Point(8) = {0.363350,0.000000,0.008866};
Point(9) = {0.335400,0.000000,0.009743};
Point(10) = {0.307450,0.000000,0.010442};
Point(11) = {0.279500,0.000000,0.010923};
Point(12) = {0.251550,0.000000,0.011158};
Point(13) = {0.223600,0.000000,0.011163};
Point(14) = {0.195650,0.000000,0.010968};
Point(15) = {0.167700,0.000000,0.010576,msTeRt};
Point(16) = {0.139750,0.000000,0.009995};
Point(17) = {0.111800,0.000000,0.009196};
Point(18) = {0.083850,0.000000,0.008156};
Point(19) = {0.055900,0.000000,0.006781};
Point(20) = {0.041925,0.000000,0.005920};
Point(21) = {0.027950,0.000000,0.004891};
Point(22) = {0.013975,0.000000,0.003617};
Point(23) = {0.006988,0.000000,0.002622};
Point(24) = {0.004193,0.000000,0.002057};
Point(25) = {0.002795,0.000000,0.001699};
Point(26) = {0.000000,0.000000,0.000000,msLeRt};
Point(27) = {0.002795,0.000000,-0.001699};
Point(28) = {0.004193,0.000000,-0.002057};
Point(29) = {0.006988,0.000000,-0.002622};
Point(30) = {0.013975,0.000000,-0.003617};
Point(31) = {0.027950,0.000000,-0.004891};
Point(32) = {0.041925,0.000000,-0.005920};
Point(33) = {0.055900,0.000000,-0.006781};
Point(34) = {0.083850,0.000000,-0.008156};
Point(35) = {0.111800,0.000000,-0.009196};
Point(36) = {0.139750,0.000000,-0.009995};
Point(37) = {0.167700,0.000000,-0.010576,msTeRt};
Point(38) = {0.195650,0.000000,-0.010968};
Point(39) = {0.223600,0.000000,-0.011163};
Point(40) = {0.251550,0.000000,-0.011158};
Point(41) = {0.279500,0.000000,-0.010923};
Point(42) = {0.307450,0.000000,-0.010442};
Point(43) = {0.335400,0.000000,-0.009743};
Point(44) = {0.363350,0.000000,-0.008866};
Point(45) = {0.391300,0.000000,-0.007837};
Point(46) = {0.419250,0.000000,-0.006680};
Point(47) = {0.447200,0.000000,-0.005406};
Point(48) = {0.475150,0.000000,-0.004075};
Point(49) = {0.503100,0.000000,-0.002739,msTeRt};
Point(50) = {0.531050,0.000000,-0.001398,msTeRt};

// Airfoil 2: agard445, 51 points 
Point(52) = {1.178128,0.762000,0.000000,msTeTp};
Point(53) = {1.159709,0.762000,0.000921};
Point(54) = {1.141290,0.762000,0.001805,msTeTp};
Point(55) = {1.122870,0.762000,0.002685};
Point(56) = {1.104451,0.762000,0.003562};
Point(57) = {1.086032,0.762000,0.004402};
Point(58) = {1.067613,0.762000,0.005165};
Point(59) = {1.049194,0.762000,0.005843};
Point(60) = {1.030775,0.762000,0.006421};
Point(61) = {1.012356,0.762000,0.006881};
Point(62) = {0.993937,0.762000,0.007198};
Point(63) = {0.975518,0.762000,0.007353};
Point(64) = {0.957099,0.762000,0.007357};
Point(65) = {0.938680,0.762000,0.007228};
Point(66) = {0.920261,0.762000,0.006970,msTeTp};
Point(67) = {0.901842,0.762000,0.006587};
Point(68) = {0.883423,0.762000,0.006060};
Point(69) = {0.865004,0.762000,0.005375};
Point(70) = {0.846585,0.762000,0.004468};
Point(71) = {0.837375,0.762000,0.003901};
Point(72) = {0.828166,0.762000,0.003223};
Point(73) = {0.818956,0.762000,0.002383};
Point(74) = {0.814351,0.762000,0.001728};
Point(75) = {0.812509,0.762000,0.001356};
Point(76) = {0.811589,0.762000,0.001120};
Point(77) = {0.809747,0.762000,0.000000,msLeTp};
Point(78) = {0.811589,0.762000,-0.001120};
Point(79) = {0.812509,0.762000,-0.001356};
Point(80) = {0.814351,0.762000,-0.001728};
Point(81) = {0.818956,0.762000,-0.002383};
Point(82) = {0.828166,0.762000,-0.003223};
Point(83) = {0.837375,0.762000,-0.003901};
Point(84) = {0.846585,0.762000,-0.004468};
Point(85) = {0.865004,0.762000,-0.005375};
Point(86) = {0.883423,0.762000,-0.006060};
Point(87) = {0.901842,0.762000,-0.006587};
Point(88) = {0.920261,0.762000,-0.006970,msTeTp};
Point(89) = {0.938680,0.762000,-0.007228};
Point(90) = {0.957099,0.762000,-0.007357};
Point(91) = {0.975518,0.762000,-0.007353};
Point(92) = {0.993937,0.762000,-0.007198};
Point(93) = {1.012356,0.762000,-0.006881};
Point(94) = {1.030775,0.762000,-0.006421};
Point(95) = {1.049194,0.762000,-0.005843};
Point(96) = {1.067613,0.762000,-0.005165};
Point(97) = {1.086032,0.762000,-0.004402};
Point(98) = {1.104451,0.762000,-0.003562};
Point(99) = {1.122870,0.762000,-0.002685};
Point(100) = {1.141290,0.762000,-0.001805,msTeTp};
Point(101) = {1.159709,0.762000,-0.000921,msTeTp};

// Box:
Point(5001) = {xO,0.000000,zO,msF};
Point(5002) = {xO+xL,0.000000,zO,msF};
Point(5003) = {xO,0.000000,zO+zL,msF};
Point(5004) = {xO+xL,0.000000,zO+zL,msF};
Point(5005) = {xO,yL,zO,msF};
Point(5006) = {xO+xL,yL,zO,msF};
Point(5007) = {xO,yL,zO+zL,msF};
Point(5008) = {xO+xL,yL,zO+zL,msF};

// Tip:
Point(5101) = {1.159709,0.762000,0.000000};
Point(5102) = {1.141290,0.762268,0.000000,2*msLeTp};
Point(5103) = {1.122870,0.762536,0.000000};
Point(5104) = {1.104451,0.762804,0.000000};
Point(5105) = {1.086032,0.763072,0.000000};
Point(5106) = {1.067613,0.763340,0.000000};
Point(5107) = {1.049194,0.763608,0.000000};
Point(5108) = {1.030775,0.763876,0.000000};
Point(5109) = {1.012356,0.764145,0.000000};
Point(5110) = {0.993937,0.764413,0.000000};
Point(5111) = {0.975518,0.764681,0.000000};
Point(5112) = {0.957099,0.764949,0.000000};
Point(5113) = {0.938680,0.765217,0.000000};
Point(5114) = {0.920261,0.765485,0.000000,2*msTeTp};
Point(5115) = {0.815077,0.767016,0.000000};

// Dummy tip center:
Point(5349) = {1.141290,0.762000,0.000000};
Point(5350) = {0.920261,0.762000,0.000000};


// Midplane:
Point(5351) = {xO+xL,0.000000,0.000000,msF};
Point(5352) = {xO+xL,0.762000,0.000000,msF};
Point(5353) = {xO+xL,yL,0.000000,msF};
Point(5354) = {1.178128,yL,0.000000,msF};
Point(5355) = {1.141290,yL,0.000000,msF};
Point(5356) = {0.920261,yL,0.000000,msF};
Point(5357) = {0.809747,yL,0.000000,msF};
Point(5358) = {xO,yL,0.000000,msF};
Point(5359) = {xO,0.762000,0.000000,msF};
Point(5360) = {xO,0.000000,0.000000,msF};

/// Lines

// Airfoil 1:
Spline(1) = {1,2,3};
Spline(2) = {3,4,5,6,7,8,9,10,11,12,13,14,15};
Spline(3) = {15,16,17,18,19,20,21,22,23,24,25,26};
Spline(4) = {26,27,28,29,30,31,32,33,34,35,36,37};
Spline(5) = {37,38,39,40,41,42,43,44,45,46,47,48,49};
Spline(6) = {49,50,1};

// Airfoil 2:
Spline(7) = {52,53,54};
Spline(8) = {54,55,56,57,58,59,60,61,62,63,64,65,66};
Spline(9) = {66,67,68,69,70,71,72,73,74,75,76,77};
Spline(10) = {77,78,79,80,81,82,83,84,85,86,87,88};
Spline(11) = {88,89,90,91,92,93,94,95,96,97,98,99,100};
Spline(12) = {100,101,52};


// Box:
Line(41) = {5001,5002};
Line(42) = {5003,5004};
Line(43) = {5005,5006};
Line(44) = {5007,5008};
Line(45) = {5001,5005};
Line(46) = {5002,5006};
Line(47) = {5003,5007};
Line(48) = {5004,5008};
Line(49) = {5001,5360};
Line(50) = {5002,5351};
Line(51) = {5003,5360};
Line(52) = {5004,5351};
Line(53) = {5005,5358};
Line(54) = {5006,5353};
Line(55) = {5007,5358};
Line(56) = {5008,5353};

// Airfoil 1 to airfoil 2:
Line(61) = {1,52};
Line(62) = {3,54};
Line(63) = {15,66};
Line(64) = {26,77};
Line(65) = {37,88};
Line(66) = {49,100};


// Tip:
Spline(121) = {52,5101,5102};
Spline(122) = {5102,5103,5104,5105,5106,5107,5108,5109,5110,5111,5112,5113,5114};
If(GMSH_MAJOR_VERSION >= 4)
    Bezier(123) = {5114,5115,77};
Else
    BSpline(123) = {5114,5115,77};
EndIf
Ellipse(124) = {54,5349,54,5102};
Ellipse(125) = {66,5350,66,5114};
Ellipse(126) = {88,5350,88,5114};
Ellipse(127) = {100,5349,100,5102};

// Midplane:
Line(131) = {5351,5352};
Line(132) = {5352,5353};
Line(133) = {5353,5354};
Line(134) = {5354,5355};
Line(135) = {5355,5356};
Line(136) = {5356,5357};
Line(137) = {5357,5358};
Line(138) = {5358,5359};
Line(139) = {5359,5360};

// Wing to midplane:
Line(161) = {1,5351};
Line(162) = {52,5352};
Line(163) = {52,5354};
Line(164) = {5102,5355};
Line(165) = {5114,5356};
Line(166) = {77,5357};
Line(167) = {77,5359};
Line(168) = {26,5360};

/// Line loops & Surfaces

// Box:
Line Loop(1) = {1,2,3,168,-51,42,52,-161};
Line Loop(2) = {48,56,-132,-131,-52};
Line Loop(3) = {44,56,133,134,135,136,137,-55};
Line Loop(4) = {47,55,138,139,-51};
Line Loop(5) = {42,48,-44,-47};
Line Loop(6) = {-6,-5,-4,168,-49,41,50,-161};
Line Loop(7) = {46,54,-132,-131,-50};
Line Loop(8) = {43,54,133,134,135,136,137,-53};
Line Loop(9) = {45,53,138,139,-49};
Line Loop(10) = {41,46,-43,-45};
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {-3};
Plane Surface(4) = {-4};
Plane Surface(5) = {-5};
Plane Surface(6) = {-6};
Plane Surface(7) = {-7};
Plane Surface(8) = {8};
Plane Surface(9) = {9};
Plane Surface(10) = {10};

// Wing 1:
Line Loop(11) = {1,62,-7,-61};
Line Loop(12) = {2,63,-8,-62};
Line Loop(13) = {3,64,-9,-63};
Line Loop(14) = {4,65,-10,-64};
Line Loop(15) = {5,66,-11,-65};
Line Loop(16) = {6,61,-12,-66};
Surface(11) = {-11};
Surface(12) = {-12};
Surface(13) = {-13};
Surface(14) = {-14};
Surface(15) = {-15};
Surface(16) = {-16};

// Wingtip:
Line Loop(71) = {7,124,-121};
Line Loop(72) = {8,125,-122,-124};
Line Loop(73) = {9,-123,-125};
Line Loop(74) = {10,126,123};
Line Loop(75) = {11,127,122,-126};
Line Loop(76) = {12,121,-127};
Surface(71) = {-71};
Surface(72) = {-72};
Surface(73) = {-73};
Surface(74) = {-74};
Surface(75) = {-75};
Surface(76) = {-76};

// Midplane:
Line Loop(81) = {161,131,-162,-61};
Line Loop(82) = {162,132,133,-163};
Line Loop(83) = {163,134,-164,-121};
Line Loop(84) = {164,135,-165,-122};
Line Loop(85) = {165,136,-166,-123};
Line Loop(86) = {167,-138,-137,-166};
Line Loop(87) = {167,139,-168,64};
Surface(81) = {81};
Surface(82) = {82};
Surface(83) = {83};
Surface(84) = {84};
Surface(85) = {85};
Surface(86) = {-86};
Surface(87) = {87};

/// Surface loops & Volumes

// Upper:
Surface Loop(1) = {11,12,13,71,72,73,1,2,3,4,5,81,82,83,84,85,86,87};
Volume(1) = {1};
// Lower:
Surface Loop(2) = {14,15,16,74,75,76,6,7,8,9,10,81,82,83,84,85,86,87};
Volume(2) = {2};



//// MESHING ALGORITHM


/// 2D:

///Wing, farfield and symmetry plane:
MeshAlgorithm Surface{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,81,82,83,84,85,86,87} = 5;

///Tip:
MeshAlgorithm Surface{71,72,73,74,75,76} = 1;

/// 3D:

Mesh.Algorithm3D = 2;
Mesh.OptimizeNetgen = 1;
Mesh.Smoothing = 10;
Mesh.SmoothNormals = 1;



//// PHYSICAL GROUPS

// Trailing edge and wake tip
Physical Line("wakeTip") = {162};
Physical Line("teTip") = {61, 162};

// Internal Field:
Physical Volume("field") = {1};
Physical Volume("field_") = {2};

// Wing:
Physical Surface("wing") = {11,12,13,71,72,73};
Physical Surface("wing_") = {14,15,16,74,75,76};

// Tip
//Physical Surface("tip") = {71,72,73,74,75,76};

// Symmetry:
Physical Surface("symmetry") = {1};
Physical Surface("symmetry_") = {6};

// Farfield:
Physical Surface("upstream") = {10};
Physical Surface("farfield") = {3,4,5,8,9};

// Downstream:
Physical Surface("downstream") = {2};
Physical Surface("downstream_") = {7};

// Wake:
Physical Surface("wake") = {81};

Coherence;
