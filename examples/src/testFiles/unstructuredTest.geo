// Gmsh project created on Thu Jan 27 13:28:48 2022
//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 10, 1, 0};
//+
MeshSize {4, 3, 2, 1} = 0.5;


//+
Curve Loop(2) = {3, 4, 1, 2};
//+
Plane Surface(2) = {2};
//+
Transfinite Surface {1} = {4, 3, 2, 1};
//+
Transfinite Curve {1, 3} = 4 Using Progression 1;
//+
Transfinite Curve {2, 4} = 4 Using Progression 1;
//+
Recombine Surface {1};
