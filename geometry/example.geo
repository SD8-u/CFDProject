//Specify mesh generation properties
Mesh.ElementOrder = 2;  //2D Quadratic Triangles
Mesh.Algorithm = 5; //Delaunay Triangulation

//Define simple square boundary
Point(1) = {0, 0, 0, 1.0};
Point(2) = {2, 0, 0, 1.0};
Point(3) = {2, 2, 0, 1.0};
Point(4) = {0, 2, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

//Tag boundary
Physical Line("Boundary") = {1, 2, 3, 4};

// Tag fluid domain insisde of square
Physical Surface("FluidDomain") = {6};