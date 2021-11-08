Mesh.MshFileVersion = 2.2;

L = 8.0;
H = 2.0;
nx = 32;
dx = H/nx;

Point(1) = {-L/2,-H/2,0.0,dx};
Point(2) = { 0.0,-H/2,0.0,dx};
Point(3) = { L/2,-H/2,0.0,dx};
Point(4) = { L/2, H/2,0.0,dx};
Point(5) = { 0.0, H/2,0.0,dx};
Point(6) = {-L/2, H/2,0.0,dx};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
Line(7) = {2, 5};

Curve Loop(1) = {1, 7, 5, 6};
Plane Surface(1) = {1};
Curve Loop(2) = {2, 3, 4, -7};
Plane Surface(2) = {2};

Transfinite Surface {1};
Transfinite Surface {2};
Recombine Surface {1, 2};

Extrude {0, 0, 2} {
  Surface{1}; Surface{2};
  Layers{32};
  Recombine;
}

Physical Surface("wall_y") = {24, 46, 16, 38};
Physical Surface("wall_z") = {1, 29, 2, 51};
Physical Volume("fluid") = {1, 2};
