//Specify mesh generation properties
Mesh.Algorithm = 5; //Delaunay Triangulation

//Define simple square boundary
Point(1) = {0, 0, 0, 1};
Point(2) = {1, 0, 0, 1};
Point(3) = {1, 1, 0, 1};
Point(4) = {0, 1, 0, 1};

Line(5) = {1, 2};
Line(6) = {2, 3};
Line(7) = {3, 4};
Line(8) = {4, 1};

Line Loop(9) = {5, 6, 7, 8};
Plane Surface(10) = {9};

//Tag boundary
Physical Line("Boundary") = {6, 7, 8};

//Tag inlet
Physical Line("Inlet") = {5};

// Tag fluid domain inside of square
Physical Surface("FluidDomain") = {10};