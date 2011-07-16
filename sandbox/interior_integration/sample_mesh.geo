cl1 = 0.25;
cl2 = 0.125;
Point(1) = {0, 0, 0, cl1};
Point(2) = {1, 0, 0, cl1};
Point(3) = {1, 1, 0, cl1};
Point(4) = {0, 1, 0, cl1};
Point(5) = {0.25, 0.25, 0, cl2};
Point(6) = {0.75, 0.25, 0, cl2};
Point(7) = {0.75, 0.75, 0, cl2};
Point(8) = {0.25, 0.75, 0, cl2};
Line(1) = {4, 3};
Line(2) = {3, 2};
Line(3) = {2, 1};
Line(4) = {1, 4};
Line(5) = {5, 8};
Line(6) = {8, 7};
Line(7) = {7, 6};
Line(8) = {6, 5};
Line Loop(11) = {6, 7, 8, 5};
Plane Surface(11) = {11};
Line Loop(12) = {1, 2, 3, 4, -5, -8, -7, -6};
Plane Surface(12) = {12};
PEC = 99;
Physical Line(PEC) = {1, 2, 3, 4};

Physical Line(100) = {5};
Physical Line(101) = {6};
Physical Line(102) = {7};
Physical Line(103) = {8};

Physical Surface(0) = {11, 12};

