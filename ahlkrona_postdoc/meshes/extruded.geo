// Unstructured base mesh extruded in the vertical (tetrahedrons)
nx = 20;
nz = 20;
h = 1000.;
l = 80000/nx;
Point(1) = {0, 0, 0, l};
Point(2) = {80000, 0, 0, l};
Point(3) = {0, 80000, 0, l};
Point(4) = {80000, 80000, 0, l};
Line(1) = {1, 2};
Line(2) = {3, 4};
Line(3) = {1, 3};
Line(4) = {2, 4};
Line Loop(5) = {1, 4, -2, -3};
Plane Surface(5) = {5};
Periodic Curve {1} = {-2};
Periodic Curve {3} = {-4};
Extrude {0,0,h} {Surface{5}; Layers{nz};}
