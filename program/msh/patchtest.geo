// ---------------------------------------------------------------------------
// patchtest.geo
//
// Unit-square patch test for the Kirchhoff plate accuracy study.
//
// Domain      : [0, a] x [0, a],  a = 1
// Approximation: structured triangular mesh (Tri3)
// Physical groups
//   Ω      : the whole surface
//   Γᵍ     : all four edges (Dirichlet boundary used by the patch test)
//
// Note: the four edges intentionally belong to one physical group only.
// Gmsh writes an element once per physical group, so putting the same edge in
// both Γᵍ and Γ¹..Γ⁴ would make ApproxOperator read duplicate boundary
// elements.
//
// Mesh density is controlled by n:
//   n = 10  ->  121 nodes, 200 triangular cells
//   n = 20  ->  441 nodes, 800 triangular cells
//   n = 40  -> 1681 nodes, 3200 triangular cells
//
// Generate the mesh with:
//   gmsh -2 program/msh/patchtest.geo -o program/msh/patchtest_10.msh
//
// To obtain an unstructured mesh instead of the transfinite grid, comment out
// the two `Transfinite ...` lines below.  The point spacing then follows
// lc = a/n.
// ---------------------------------------------------------------------------

a  = 1.0;
n  = 64;
lc = a/n;

Point(1) = {0.0, 0.0, 0.0, lc};
Point(2) = {  a, 0.0, 0.0, lc};
Point(3) = {  a,   a, 0.0, lc};
Point(4) = {0.0,   a, 0.0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Structured grid: n intervals per edge, (n+1)^2 nodes.
Transfinite Curve{1, 2, 3, 4} = n + 1;
Transfinite Surface{1};

// Full boundary: the Dirichlet group used by the patch test.
Physical Curve("Γᵍ") = {1, 2, 3, 4};

Physical Surface("Ω") = {1};

Mesh.MshFileVersion      = 2;
Mesh.ElementOrder        = 1;
Mesh.Renumber            = 0;
Mesh.Algorithm           = 6;
Mesh 2;
