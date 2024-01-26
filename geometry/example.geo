//Specify mesh generation properties
Mesh.Algorithm = 5; //Delaunay Triangulation

//Define simple square boundary
Point(1) = {0, 0, 0, 1.0};
Point(2) = {10, 0, 0, 1.0};
Point(3) = {10, 10, 0, 1.0};
Point(4) = {0, 10, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

//Tag boundary
Physical Line("Boundary") = {2, 3, 4};

//Tag inlet
Physical Line("Inlet") = {1};

// Tag fluid domain insisde of square
Physical Surface("FluidDomain") = {6};