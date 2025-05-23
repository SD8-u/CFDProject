// Backwards facing step geometry
//+
Point(1) = {-0.5, 0, 0, 1.0};
//+
Point(2) = {0.5, 0, 0, 4.0};
//+
Point(3) = {0.5, -0.5, -0, 1.0};

Point(4) = {5, -0.5, -0, 1.0};

Point(5) = {5, 0.5, -0, 1.0};

Point(6) = {-0.5, 0.5, -0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

Line Loop(1) = {1, 2, 3, 4, 5, 6};

Plane Surface(1) = {1};

//Tag boundary
Physical Line("Boundary") = {1, 2, 3, 5};

//Tag inlet
Physical Line("Inlet") = {6};

//Tag outlet
Physical Line("Outlet") = {4};

Physical Surface("FluidDomain") = {1};