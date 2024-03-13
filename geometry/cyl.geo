// Gmsh project created on Wed Mar 13 18:40:41 2024
//+
Point(1) = {-1, 0, 0, 1.0};
//+
Point(2) = {-2, -1.5, 0, 1.0};
//+
Point(3) = {-2, 1.5, 0, 1.0};
//+
Point(4) = {5, -1.5, 0, 1.0};
//+
Point(5) = {5, 1.5, 0, 1.0};
//+
Point(6) = {-1, 0.9, 0, 1.0};
//+
Point(7) = {-1, -0.9, 0, 1.0};
//+
Circle(1) = {6, 1, 7};
//+
Circle(2) = {7, 1, 6};
//+
Line(3) = {5, 3};
//+
Line(4) = {2, 3};
//+
Line(5) = {4, 2};
//+
Line(6) = {5, 4};
//+
Curve Loop(1) = {3, -6, -5, -4};
//+
Curve Loop(2) = {2, 1};
//+
Plane Surface(1) = {1, 2};

//Tag boundary
Physical Line("Boundary") = {1, 2, 3, 5};

//Tag inlet
Physical Line("Inlet") = {4};

//Tag outlet
Physical Line("Outlet") = {6};

// Tag fluid domain inside of square
Physical Surface("FluidDomain") = {1};