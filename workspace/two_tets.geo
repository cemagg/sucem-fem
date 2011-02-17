h=1;
a=0.5;
Point(1) = {-a, a,-a,h};
Point(2) = { a, a, a,h};
Point(3) = { a,-a,-a,h};
Point(4) = {-a,-a, a,h};
Point(5) = {-a,-a,-a,h};
Line(1) = {1,2};
Line(2) = {1,3};
Line(3) = {1,4};
Line(4) = {2,3};
Line(5) = {2,4};
Line(6) = {3,4};
Line(7) = {1,5};
Line(8) = {3,5};
Line(9) = {4,5};
Line Loop(10) = {1,4,-2};
Plane Surface(11) = {10};
Line Loop(12) = {1,5,-3};
Plane Surface(13) = {12};
Line Loop(14) = {2,6,-3};
Plane Surface(15) = {14};
Line Loop(16) = {4,6,-5};
Plane Surface(17) = {16};
Line Loop(18) = {2,8,-7};
Plane Surface(19) = {18};
Line Loop(20) = {3,9,-7};
Plane Surface(21) = {20};
Line Loop(22) = {6,9,-8};
Plane Surface(23) = {22};
Surface Loop(24) = {11,-13,-17,15};
Volume(25) = {24};
Surface Loop(26) = {15,-19,-21,23};
Volume(27) = {26};

fc1 = newp; Point(fc1) = {0.16666666666667,0.16666666666667,-0.16666666666667,h};
np1 = newp; Point(np1) = {0.74401693585629,0.74401693585629,-0.74401693585629,h};
n1 = newl; Line(n1) = {fc1, np1};

fc2 = newp; Point(fc2) = {-0.16666666666667,0.16666666666667,0.16666666666667,h};
np2 = newp; Point(np2) = {-0.74401693585629,0.74401693585629,0.74401693585629,h};
n2 = newl; Line(n2) = {fc2, np2};

fc3 = newp; Point(fc3) = {0.16666666666667,-0.16666666666667,0.16666666666667,h};
np3 = newp; Point(np3) = {0.74401693585629,-0.74401693585629,0.74401693585629,h};
n3 = newl; Line(n3) = {fc3, np3};

fc4 = newp; Point(fc4) = {-0.16666666666667,-0.16666666666667,-0.5,h};
np4 = newp; Point(np4) = {-0.16666666666667,-0.16666666666667,-1.5,h};
n4 = newl; Line(n4) = {fc4, np4};

fc5 = newp; Point(fc5) = {-0.5,-0.16666666666667,-0.16666666666667,h};
np5 = newp; Point(np5) = {-1.5,-0.16666666666667,-0.16666666666667,h};
n5 = newl; Line(n5) = {fc5, np5};

fc6 = newp; Point(fc6) = {-0.16666666666667,-0.5,-0.16666666666667,h};
np6 = newp; Point(np6) = {-0.16666666666667,-1.5,-0.16666666666667,h};
n6 = newl; Line(n6) = {fc6, np6};


