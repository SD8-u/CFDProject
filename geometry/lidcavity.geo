//Specify mesh generation properties
Mesh.Algorithm = 5; //Delaunay Triangulation

//Define simple square boundary
Point(1) = {0, 0, 0, 2};
Point(2) = {1, 0, 0, 2};
Point(3) = {1, 1, 0, 2};
Point(4) = {0, 1, 0, 2};

Line(1) = {2, 3};
Line(2) = {3, 4};
Line(3) = {4, 1};
Line(4) = {1, 2};

Line Loop(5) = {2, 3, 4, 1};
Plane Surface(6) = {5};

//Tag boundary
Physical Line("Boundary") = {1, 4, 3};

//Tag inlet
Physical Line("Inlet") = {2};

// Tag fluid domain inside of square
Physical Surface("FluidDomain") = {6};