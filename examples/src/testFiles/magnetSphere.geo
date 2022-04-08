//+

//Ellipse(1) = {0, 0, 0, 0.5, 0.25, 0, 2*Pi};
//+
//Ellipse(2) = {0, 0, 0, 3, 3, 0, 2*Pi};
//+
//MeshSize {1} = 0.1;
//+
//+
//Curve Loop(1) = {2};
//+
//Plane Surface(1) = {1};


SetFactory("OpenCASCADE");
hv = 1.0;
Sphere(1) = {0.0,0.0,hv/2, 3.0, -Pi/2, Pi/2, 2*Pi};

Disk(2) = {0, 0, 0, 0.5, 0.5};



Printf("%f",hv);
zyl[] = Extrude {0, 0, hv} {
  Surface{2};
};

//+
MeshSize {3, 4} = 0.1;
MeshSize {1, 2} = 1;


//Surface{2} In Volume{1};
//Surface{3} In Volume{1};
//Surface{4} In Volume{1};
//Surface{2, 3, 1} In Volume{2};

//Surface Loop(301) = {2,3,4};
//Volume(302) = {301};

//resulting_vol[] = BooleanUnion{ Volume{302}; Delete;}{ Volume{1}; Delete;};

BooleanFragments{ Volume{1}; }{Surface{2,3,4};}
//+
Recursive Delete {
  Volume{1};
  Volume{4};
}
