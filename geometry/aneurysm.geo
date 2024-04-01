//Anyeursym geometry
//Specify mesh generation properties
Mesh.Algorithm = 5; //Delaunay Triangulation

Point(1) = {0, 0, 0, 1.0};
Point(2) = {3, 0, 0, 1.0};
Point(3) = {3, 0.25, 0, 1.0};
Point(4) = {0, 0.25, 0, 1.0};

Point(7) = {2.5, 0, 0, 1.0};
Point(8) = {2, 0, 0, 1.0};
Point(9) = {1.5, 0, 0, 1.0};
Point(10) = {2, -0.5, 0, 1.0};

Line(1) = {1, 9};
Circle(11) = {9, 8, 10};
Circle(12) = {10, 8, 7};
Line(13) = {7, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(5) = {1, 11, 12, 13, 2, 3, 4};
Plane Surface(6) = {5};

//Tag boundary
Physical Line("Boundary") = {1, 11, 12, 13, 3};

//Tag inlet
Physical Line("Inlet") = {2};

//Tag outlet
Physical Line("Outlet") = {4};

// Tag fluid domain inside of square
Physical Surface("FluidDomain") = {6};