// Unstructured base mesh extruded in the vertical (tetrahedrons)
// Normally, the vertical resolution is much higher than horizontal;
// hence the domain must be scaled in the vertical __before__ meshing.
// Typically, the resolution is around 50 m in the vertical and ~2km in the horizontal.
// For the ISMIP HOM A domain (80x80x1 km), this would give a vertical scaling factor of 40.

nz = 20; // approximate number of cells in vertical
nx = 40; // approximate number of cells in horizontal
h = 1000.; // physical height of domain (what it should be in the simulation)
L = 80000; // physical length of domain
lc = L / nx; // characteristic length
scaling_factor = lc *nz / h;
H = h * scaling_factor;
Point(1) = {0, 0, 0, lc};
Point(2) = {L, 0, 0, lc};
Point(3) = {0, L, 0, lc};
Point(4) = {L, L, 0, lc};
Point(1) = {0, 0, 0, lc};
Point(2) = {L, 0, 0, lc};
Point(3) = {0, L, 0, lc};
Point(4) = {L, L, 0, lc};

Line(1) = {1, 2};
Line(2) = {3, 4};
Line(3) = {1, 3};
Line(4) = {2, 4};
Line Loop(5) = {1, 4, -2, -3};
Plane Surface(5) = {5};
Periodic Curve {1} = {-2};
Periodic Curve {3} = {-4};
extr_out[] = Extrude {0,0,H} {Surface{5}; };
// Periodic Surface {extr_out[2]} = {extr_out[4]};
// Periodic Surface {extr_out[3]} = {extr_out[5]};
Printf("extr_out[0] = %g", extr_out[0]);
Printf("extr_out[1] = %g", extr_out[1]);
Printf("extr_out[2] = %g", extr_out[2]);
Printf("extr_out[3] = %g", extr_out[3]);
Printf("extr_out[4] = %g", extr_out[4]);
Printf("extr_out[5] = %g", extr_out[5]);
