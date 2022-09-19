// Gmsh project created on Wed Jun 15 11:47:52 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 44, 0, 1.0};
//+
Point(3) = {48, 60, 0, 1.0};
//+
Point(4) = {48, 44, 0, 1.0};
//+
Point(5) = {0, 0, 1, 1.0};
//+
Point(6) = {0, 44, 1, 1.0};
//+
Point(7) = {48, 60, 1, 1.0};
//+
Point(8) = {48, 44, 1, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 5};
//+
Line(9) = {1, 5};
//+
Line(10) = {2, 6};
//+
Line(11) = {3, 7};
//+
Line(12) = {4, 8};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Curve Loop(2) = {5, 6, 7, 8};
//+
Curve Loop(3) = {9, 1, 10, 5};
//+
Curve Loop(4) = {10, 2, 11, 6};
//+
Curve Loop(5) = {11, 3, 12, 7};
//+
Curve Loop(6) = {12, 4, 9, 8};
//+
Plane Surface(1) = {1};
//+
Plane Surface(2) = {2};
//+
Plane Surface(3) = {3};
//+
Plane Surface(4) = {4};
//+
Plane Surface(5) = {5};
//+
Plane Surface(6) = {6};
//+
Surface Loop(1) = {3, 6, 5, 4, 2, 1};
//+
Volume(1) = {1};
//+
Transfinite Volume {1};
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Surface {3};
//+
Transfinite Surface {4};
//+
Transfinite Surface {5};
//+
Transfinite Surface {6};
//+
Transfinite Curve {1, 3} = 3 Using Progression 1;
//+
Transfinite Curve {2, 4} = 3 Using Progression 1;
//+
Transfinite Curve {5, 7} = 3 Using Progression 1;
//+
Transfinite Curve {6, 8} = 3 Using Progression 1;
//+
Transfinite Curve {9, 10} = 3 Using Progression 1;
//+
Transfinite Curve {10, 11} = 3 Using Progression 1;
//+
Transfinite Curve {11, 12} = 3 Using Progression 1;
//+
Transfinite Curve {12, 9} = 3 Using Progression 1;
//+
Recombine Surface {1};
//+
Recombine Surface {2};
//+
Recombine Surface {3};
//+
Recombine Surface {4};
//+
Recombine Surface {5};
//+
Recombine Surface {6};
//+
Recombine Volume {1};
