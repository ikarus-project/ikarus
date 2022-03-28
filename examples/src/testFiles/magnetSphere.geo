//+
SetFactory("OpenCASCADE");
Ellipse(1) = {0, 0, 0, 0.5, 0.25, 0, 2*Pi};
//+
Ellipse(2) = {0, 0, 0, 3, 3, 0, 2*Pi};
//+
MeshSize {1} = 0.1;
//+
//+
Curve Loop(1) = {2};
//+
Plane Surface(1) = {1};
