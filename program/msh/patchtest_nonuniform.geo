// ---------------------------------------------------------------------------
// patchtest_nonuniform.geo
//
// Non-uniform (graded, unstructured) version of patchtest.geo.
//
// Domain      : [0, a] x [0, a],  a = 1
// Approximation: unstructured triangular mesh (Tri3)
// Physical groups
//   Ω      : the whole surface
//   Γᵍ     : all four edges (Dirichlet boundary used by the patch test)
//
// The mesh size is prescribed at the four corners and interpolated by Gmsh:
//
//   (0,0) : lc1 = 0.06   (finest)
//   (1,0) : lc2 = 0.10
//   (1,1) : lc3 = 0.14   (coarsest)
//   (0,1) : lc2 = 0.10
//
// This produces a reproducible non-uniform node distribution without random
// perturbation.  The default mesh has about 169 nodes and 340 triangles.
//
// Generate the mesh with:
//   gmsh -2 program/msh/patchtest_nonuniform.geo \
//        -o program/msh/patchtest_nonuniform.msh
//
// To change the grading, edit lc1, lc2 and lc3 below.
// ---------------------------------------------------------------------------

a   = 1.0;
lc1 = 0.06;
lc2 = 0.10;
lc3 = 0.14;

Point(1) = {0.0, 0.0, 0.0, lc1};
Point(2) = {  a, 0.0, 0.0, lc2};
Point(3) = {  a,   a, 0.0, lc3};
Point(4) = {0.0,   a, 0.0, lc2};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Full boundary: the Dirichlet group used by the patch test.
Physical Curve("Γᵍ") = {1, 2, 3, 4};

Physical Surface("Ω") = {1};

Mesh.MshFileVersion      = 2;
Mesh.ElementOrder        = 1;
Mesh.Renumber            = 0;
Mesh.Algorithm           = 6;
Mesh 2;
