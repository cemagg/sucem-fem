// Gmsh project created on Mon Apr  4 13:01:19 2011
a = 1;
b = 4*a/9;
c = 5*a/9;
s = 5*a/18;
t = a/2;
u = 2*a/9;
d = (a-t)/2;
w = (c-u)/2;

lc0 = a/2;
delta_m = 0;
eps_r = 2.05;


Point(1) = { 0, c, 0, lc0};
Point(2) = { 0, c, b, lc0};
Point(3) = { 0, 0, b, lc0};
Point(4) = { 0, 0, 0, lc0};

Point(5) = { d, w+u, 0, lc0};
Point(6) = { d, w+u, s, lc0};
Point(7) = { d, w  , s, lc0};
Point(8) = { d, w  , 0, lc0};

Point(9) =  { t+d, w+u, 0, lc0};
Point(10) = { t+d, w+u, s, lc0};
Point(11) = { t+d, w  , s, lc0};
Point(12) = { t+d, w  , 0, lc0};

Point(13) = { a, c, 0, lc0};
Point(14) = { a, c, b, lc0};
Point(15) = { a, 0, b, lc0};
Point(16) = { a, 0, 0, lc0};

// Mesh control
Field[1] = Box;
Field[1].VIn = lc0 / Sqrt(eps_r);
Field[1].VOut = lc0;
Field[1].XMin = d - delta_m; 
Field[1].XMax = d+t + delta_m;
Field[1].YMin = w - delta_m;
Field[1].YMax = w+u + delta_m;
Field[1].ZMin = 0;
Field[1].ZMax = s;
Background Field = 1;

Line(1) = {13, 14};
Line(2) = {14, 15};
Line(3) = {15, 16};
Line(4) = {16, 13};
Line(5) = {9, 10};
Line(6) = {10, 11};
Line(7) = {11, 12};
Line(8) = {12, 9};
Line(9) = {9, 5};
Line(10) = {5, 8};
Line(11) = {8, 12};
Line(12) = {11, 7};
Line(13) = {7, 8};
Line(14) = {7, 6};
Line(15) = {6, 10};
Line(16) = {6, 5};
Line(17) = {1, 4};
Line(18) = {4, 3};
Line(19) = {3, 2};
Line(20) = {2, 1};
Line(21) = {16, 4};
Line(22) = {3, 15};
Line(23) = {2, 14};
Line(24) = {1, 13};
Line Loop(25) = {21, -17, 24, -4};
Line Loop(26) = {9, 10, 11, 8};
Plane Surface(27) = {25, 26};
Plane Surface(28) = {26};
Line Loop(29) = {21, 18, 22, 3};
Plane Surface(30) = {29};
Line Loop(31) = {11, -7, 12, 13};
Plane Surface(32) = {31};
Line Loop(33) = {18, 19, 20, 17};
Plane Surface(34) = {33};
Line Loop(35) = {14, 16, 10, -13};
Plane Surface(36) = {35};
Line Loop(37) = {5, 6, 7, 8};
Plane Surface(38) = {37};
Line Loop(39) = {3, 4, 1, 2};
Plane Surface(40) = {39};
Line Loop(41) = {22, -2, -23, -19};
Plane Surface(42) = {41};
Line Loop(43) = {6, 12, 14, 15};
Plane Surface(44) = {43};
Line Loop(45) = {15, -5, 9, -16};
Plane Surface(46) = {45};
Line Loop(47) = {23, -1, -24, -20};
Plane Surface(48) = {47};
Surface Loop(49) = {42, 30, 27, 34, 48, 40, 44, 38, 46, 36, 32};
Volume(50) = {49};
Surface Loop(51) = {38, 46, 44, 32, 36, 28};
Volume(52) = {51};

Physical Volume(1000) = {50};
Physical Volume(1001) = {52};