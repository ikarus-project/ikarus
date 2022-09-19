// Gmsh project created on Wed Jun 15 11:47:52 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 29.0};
//+
Point(2) = {0, 44, 0, 18.0};
//+
Point(3) = {48, 60, 0, 25.0};
//+
Point(4) = {48, 44, 0, 10.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Recombine Surface {1};
