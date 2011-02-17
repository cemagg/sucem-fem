lc = 1/3;
ro = 1;
ri = 2/3;
len = 10;
Point(1) = {0,0,0,lc};
Point(2) = {0,ro,0,lc};
Point(3) = {0,-ro,0,lc};
// Point(4) = {ro,0,0,lc};
// Point(5) = {-ro,0,0,lc};
Point(6) = {0,ri,0,lc};
Point(7) = {0,-ri,0,lc};
// Point(8) = {ri,0,0,lc};
// Point(9) = {-ri,0,0,lc};
Circle(1) = {2,1,3};
Circle(2) = {6,1,7};
// Circle(2) = {3,1,5};
// Circle(3) = {5,1,2};
// Circle(4) = {2,1,4};
// Circle(5) = {8,1,7};
// Circle(6) = {7,1,9};
// Circle(7) = {9,1,6};
// Circle(8) = {6,1,8};
// Line(9) = {2,6};
// Line(10) = {4,8};
// Line(11) = {3,7};
// Line(12) = {5,9};
// Line Loop(13) = {4,10,-8,-9};
// Plane Surface(14) = {13};
// Line Loop(15) = {1,11,-5,-10};
// Plane Surface(16) = {15};
// Line Loop(17) = {11,6,-12,-2};
// Plane Surface(18) = {17};
// Line Loop(19) = {12,7,-9,-3};
// Plane Surface(20) = {19};
Line(3) = {3,7};
Line(4) = {6,2};
Line Loop(5) = {1,3,-2,4};
Plane Surface(6) = {5};
// Translate {0,0,len/2} {
//   Duplicata { Surface{20,14,16,18}; }
// }
Translate {0,0,len/2} {
  Duplicata { Surface{6}; }
}
// Translate {0,0,len} {
//   Duplicata { Surface{20,14,16,18}; }
// }
Translate {0,0,len} {
  Duplicata { Surface{6}; }
}
// Line(41) = {6,16};
// Line(42) = {20,2};
// Line(43) = {34,8};
// Line(44) = {4,30};
// Line(45) = {52,7};
// Line(46) = {3,48};
// Line(47) = {11,9};
// Line(48) = {5,10};
Mesh.Algorithm3D = 5; // Tetgen 79k els, rho:  0.0432302595079 1/rho: 23.1319453407
//Mesh.Algorithm3D = 4; // netgen 50k els, rho:  0.0425506141222 1/rho: 23.5014234372
//Mesh.Algorithm3D = 1; // delaunay 79k els, rho:  0.0189926852938 1/rho: 52.6518490951

Line(17) = {20,8};
Line(18) = {31,19};
Line(19) = {26,14};
Line(20) = {22,10};
Line(21) = {8,2};
Line(22) = {19,6};
Line(23) = {14,7};
Line(24) = {10,3};
Line Loop(25) = {17,-11,-18,16};
Plane Surface(26) = {25};
Line Loop(27) = {21,-4,-22,11};
Plane Surface(28) = {27};
Line Loop(29) = {19,-9,-20,14};
Plane Surface(30) = {29};
Line Loop(31) = {23,-3,-24,9};
Plane Surface(32) = {31};
Line Loop(33) = {8,-20,-13,17};
Ruled Surface(34) = {33};
Line Loop(35) = {18,-10,-19,15};
Ruled Surface(36) = {35};
Line Loop(37) = {22,2,-23,10};
Ruled Surface(38) = {37};
Line Loop(39) = {1,-24,-8,21};
Ruled Surface(40) = {39};
Surface Loop(41) = {40,6,32,38,28,7};
Volume(42) = {41};
Surface Loop(43) = {26,34,30,36,12,7};
Volume(44) = {43};
