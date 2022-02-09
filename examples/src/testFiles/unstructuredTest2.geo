// Gmsh project created on Fri Jan 28 09:00:37 2022
//SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 1, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {1, 0, 0, 1.0};
//+
Line(1) = {2, 3};
//+
Line(2) = {3, 4};
//+
Line(3) = {4, 1};
//+
Line(4) = {1, 2};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
//Transfinite Surface {1} = {2, 3, 4, 1};
//+
//Transfinite Curve {4, 2} = 3 Using Progression 1;
//+
//Transfinite Curve {1, 3} = 3 Using Progression 1;
//+
MeshSize {2, 3, 4, 1} = 0.1;
