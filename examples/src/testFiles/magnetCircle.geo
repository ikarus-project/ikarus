// Gmsh project created on Fri Mar 25 18:29:31 2022
SetFactory("OpenCASCADE");
//+
Circle(1) = {-0, -0, -0, 0.5, 0, 2*Pi};
//+
Disk(2) = {-0, -0, 0, 3, 3};
//+
MeshSize {1} = 0.1;
//+
Curve{1} In Surface{2};
//+

