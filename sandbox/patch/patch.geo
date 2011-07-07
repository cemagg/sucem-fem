Lx = 31.18e-3;
Lx = 31.18e-3;
Ly = 46.75e-3;
Sx = Lx*2;
Sy = Ly*2;
Dp = 1.3e-3;
IIf = 8.9e-3;
Hs = 2.87e-3;
eps_r = 2.2;
freq = 3e9;
lam = 299792458/freq;
Lc = lam/8;

Bx = Sx+lam;
By = Sy+lam;
Bz = Hs+lam;

// Patch Points
Point(1) = {-Lx/2, -Ly/2, Hs, Lc};
Point(2) = {-Lx/2, Ly/2, Hs, Lc};
Point(3) = {Lx/2, Ly/2, Hs, Lc};
Point(4) = {Lx/2, -Ly/2, Hs, Lc};
// Substrate Points
Point(5) = {-Sx/2, -Sy/2, 0, Lc};
Point(6) = {-Sx/2, Sy/2, 0, Lc};
Point(7) = {Sx/2, Sy/2, 0, Lc};
Point(8) = {Sx/2, -Sy/2, 0, Lc};
Point(9) = {-Sx/2, -Sy/2, Hs, Lc};
Point(10) = {-Sx/2, Sy/2, Hs, Lc};
Point(11) = {Sx/2, Sy/2, Hs, Lc};
Point(12) = {Sx/2, -Sy/2, Hs, Lc};
// Surrounding Domain Points
Point(13) = {-Bx/2, -By/2, -Bz/2+Hs/2, Lc};
Point(14) = {-Bx/2, By/2, -Bz/2+Hs/2, Lc};
Point(15) = {Bx/2, By/2, -Bz/2+Hs/2, Lc};
Point(16) = {Bx/2, -By/2, -Bz/2+Hs/2, Lc};
Point(17) = {-Bx/2, -By/2, Bz/2+Hs/2, Lc};
Point(18) = {-Bx/2, By/2, Bz/2+Hs/2, Lc};
Point(19) = {Bx/2, By/2, Bz/2+Hs/2, Lc};
Point(20) = {Bx/2, -By/2, Bz/2+Hs/2, Lc};
Line(1) = {1, 4};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 1};
Line(5) = {1, 9};
Line(6) = {9, 12};
Line(7) = {12, 4};
Line(8) = {12, 11};
Line(9) = {11, 3};
Line(10) = {11, 10};
Line(11) = {10, 2};
Line(12) = {10, 9};
Line(13) = {10, 18};
Line(14) = {18, 17};
Line(15) = {17, 9};
Line(16) = {17, 20};
Line(17) = {20, 19};
Line(18) = {19, 18};
Line(19) = {20, 12};
Line(20) = {11, 19};
Line(21) = {9, 5};
Line(22) = {5, 6};
Line(23) = {6, 10};
Line(24) = {5, 8};
Line(25) = {8, 12};
Line(26) = {8, 7};
Line(27) = {7, 11};
Line(28) = {7, 6};
Line(29) = {7, 15};
Line(30) = {6, 14};
Line(31) = {14, 15};
Line(32) = {15, 16};
Line(33) = {16, 13};
Line(34) = {13, 14};
Line(35) = {16, 8};
Line(36) = {13, 5};
Line(37) = {17, 13};
Line(38) = {14, 18};
Line(39) = {20, 16};
Line(40) = {15, 19};
// Patch surface
Line Loop(41) = {4, 1, 2, 3};
Plane Surface(42) = {41};
// Substrate top surfaces
Line Loop(43) = {11, -3, -9, 10};
Plane Surface(44) = {43};
Line Loop(45) = {9, -2, -7, 8};
Plane Surface(46) = {45};
Line Loop(47) = {7, -1, 5, 6};
Plane Surface(48) = {47};
Line Loop(49) = {5, -12, 11, 4};
Plane Surface(50) = {49};
//Substrate side surfaces
Line Loop(51) = {12, 21, 22, 23};
Plane Surface(52) = {51};
Line Loop(53) = {21, 24, 25, -6};
Plane Surface(54) = {53};
Line Loop(55) = {25, 8, -27, -26};
Plane Surface(56) = {55};
Line Loop(57) = {23, -10, -27, 28};
Plane Surface(58) = {57};
// Substrate bottom surface
Line Loop(95) = {22, -28, -26, -24};
Plane Surface(96) = {95};
// Top diagonal air surfaces
Line Loop(59) = {12, -15, -14, -13};
Plane Surface(60) = {59};
Line Loop(61) = {15, 6, -19, -16};
Plane Surface(62) = {61};
Line Loop(63) = {19, 8, 20, -17};
Plane Surface(64) = {63};
Line Loop(65) = {13, -18, -20, 10};
Plane Surface(66) = {65};
// Bottom diagonal air surfaces
Line Loop(67) = {30, 31, -29, 28};
Plane Surface(68) = {67};
Line Loop(69) = {29, 32, 35, 26};
Plane Surface(70) = {69};
Line Loop(71) = {30, -34, 36, 22};
Plane Surface(72) = {71};
Line Loop(73) = {35, -24, -36, -33};
Plane Surface(74) = {73};
// Diagonal surfaces air between top and bottom
Line Loop(75) = {21, -36, -37, 15};
Plane Surface(76) = {75};
Line Loop(77) = {35, 25, -19, 39};
Plane Surface(78) = {77};
Line Loop(79) = {27, 20, -40, -29};
Plane Surface(80) = {79};
Line Loop(81) = {13, -38, -30, 23};
Plane Surface(82) = {81};
// Bounding Box Surfaces
Line Loop(83) = {16, 17, 18, 14};
Plane Surface(84) = {83};
Line Loop(85) = {18, -38, 31, 40};
Plane Surface(86) = {85};
Line Loop(87) = {38, 14, 37, 34};
Plane Surface(88) = {87};
Line Loop(89) = {40, -17, 39, -32};
Plane Surface(90) = {89};
Line Loop(91) = {39, 33, -37, 16};
Plane Surface(92) = {91};
Line Loop(93) = {34, 31, 32, 33};
Plane Surface(94) = {93};
// Substrate Volume
Surface Loop(97) = {56, 46, 44, 50, 48, 42, 58, 52, 54, 96};
Volume(98) = {97};
// Volume above patch
Surface Loop(99) = {42, 50, 48, 46, 44, 62, 60, 66, 64, 84};
Volume(100) = {99};
// Volume below patch 
Surface Loop(101) = {70, 68, 72, 96, 74, 94};
Volume(102) = {101};
// Volume in +ve x direction from patch
Surface Loop(103) = {90, 78, 56, 64, 80, 70};
Volume(104) = {103};
// Volume in -ve x direction from patch
Surface Loop(105) = {52, 82, 60, 76, 72, 88};
Volume(106) = {105};
// Volume in +ve y direction from patch
Surface Loop(107) = {92, 62, 76, 54, 78, 74};
Volume(108) = {107};
// Volume in -ve y direction from patch
Surface Loop(109) = {58, 80, 66, 86, 68, 82};
Volume(110) = {109};
// Patch physical metal surfaces
Physical Surface(10) = {42, 96};
// Substrate physical volume
Physical Volume(20) = {98};
// Air physical volumes
Physical Volume(0) = {110, 106, 104, 108, 102, 100};
Characteristic Length {28, 1, 2, 3, 4, 5,6, 7, 8, 9, 10, 11, 12, 21, 22, 23, 24, 25, 26, 27} = Lc/2.5/Sqrt(eps_r);