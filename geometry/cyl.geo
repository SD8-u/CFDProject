//Flow around a cylinder geometry
//+
Point(1) = {0.5, 0, 0, 1.0};
//+
Point(2) = {-0.5, -0.675, 0, 1.0};
//+
Point(3) = {-0.5, 0.675, 0, 1.0};
//+
Point(4) = {1.95, -0.675, 0, 1.0};
//+
Point(5) = {1.95, 0.675, 0, 1.0};
//+
Point(6) = {0.5, 0.225, 0, 1.0};
//+
Point(7) = {0.5, -0.225, 0, 1.0};
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