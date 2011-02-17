h=1;
a=0.5;
b=1/6;
Point(1) = {-a, a,-a,h};
Point(2) = {-a,-a, a,h};
Point(3) = {-b, b, b,h};
Point(4) = {-b,-b,-b,h};
Line(1) = {1,2};
Line(2) = {1,3};
Line(3) = {1,4};
Line(4) = {2,3};
Line(5) = {2,4};
Line(6) = {3,4};
Line Loop(10) = {1,4,-2};
Plane Surface(11) = {10};
Line Loop(12) = {1,5,-3};
Plane Surface(13) = {12};
Line Loop(14) = {2,6,-3};
Plane Surface(15) = {14};
Line Loop(16) = {4,6,-5};
Plane Surface(17) = {16};
Surface Loop(24) = {11,-13,-17,15};
Volume(25) = {24};

fc1 = newp; Point(fc1) = {(-b-2*a)/3,b/3 , b/3, h};
fc2 = newp; Point(fc2) = {(-b-2*a)/3,-b/3,-b/3, h};
fc3 = newp; Point(fc3) = {(-2*b-a)/3,a/3 ,-a/3, h};
fc4 = newp; Point(fc4) = {(-2*b-a)/3,-a/3, a/3, h};

np1 = newp; Point(np1) = {(-b-2*a)/3-1/3,  b/3+1/3,  b/3+1/3,h};
np2 = newp; Point(np2) = {(-b-2*a)/3-1/3, -b/3-1/3, -b/3-1/3,h};
np3 = newp; Point(np3) = {(-2*b-a)/3+1/3,  a/3+1/9, -a/3-1/9,h};
np4 = newp; Point(np4) = {(-2*b-a)/3+1/3, -a/3-1/9,  a/3+1/9,h};

n1 = newl; Line(n1) = {fc1, np1};
n2 = newl; Line(n2) = {fc2, np2};
n3 = newl; Line(n3) = {fc3, np3};
n4 = newl; Line(n4) = {fc4, np4};

// View "normal 1" {
//   VP (-0.5,-0.5,-0.5){0,0,1};};
