# ---------------------------------------------------------------------------
# test_rkgsi_consistency.jl
#
# Direct test of the RKGSI variational consistency of the assembled mixed
# Hellinger-Reissner system.
#
# The exact polynomial
#
#     w = (x + y + 1)^p
#
# is inserted into the discrete equations.  If the RKGSI operators are
# variationally consistent, the residual of the exact solution must vanish.
#
# The exact M coefficients are the monomial coefficients of
#
#     M11 = M22 = -D p(p-1)(1+ν) (x+y+1)^(p-2)
#     M12       = -D p(p-1)(1-ν) (x+y+1)^(p-2)
#
# The exact w coefficients are the nodal values w(x_I), because the RK space
# reproduces polynomials exactly.
#
# Run from the project root:
#
#     julia program/test_rkgsi_consistency.jl
# ---------------------------------------------------------------------------

using ApproxOperator
import ApproxOperator: AbstractElement, get𝑛𝑝
import ApproxOperator.GmshImport:
    getPhysicalGroups, get𝑿ᵢ, getElements,
    getPiecewiseElements, getPiecewiseBoundaryElements
import ApproxOperator.KirchhoffPlate:
    ∫MC⁻¹MdΩ, ∫∇²MwdΩ,
    ∫VwdΓ, ∫MθdΓ, wP,
    ∫wVdΓ, ∫θₙMₙₙdΓ, wΔMₙₛ,
    ∫wqdΩ
import Gmsh: gmsh

using LinearAlgebra
using Printf
using CairoMakie
set_theme!(fonts = (
    regular = "Times New Roman",
    bold = "Times New Roman Bold",
    italic = "Times New Roman Italic",
    bolditalic = "Times New Roman Bold Italic",
))

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
const p        = 3        # RK basis degree
const n        = 2        # variational-consistency order / polynomial degree
const q        = 2        # M basis degree: 0 constant, 1 linear, 2 quadratic
const meshfile = joinpath(@__DIR__, "msh", "patchtest_rkgsi.msh")
const D        = 2.0e4
const ν        = 0.3
const order    = p == 2 ? 4 : p == 3 ? 6 : 8
const plotfile = joinpath(@__DIR__, "figures", "rkgsi_consistency_surface.png")

# ---------------------------------------------------------------------------
# Manufactured solution
# ---------------------------------------------------------------------------
@inline g(x, y) = x + y + 1.0

w(x, y, z) = g(x, y)^n
w₁(x, y, z) = n*g(x, y)^abs(n-1)
w₂(x, y, z) = n*g(x, y)^abs(n-1)
w₁₁(x, y, z) = n*(n-1)*g(x, y)^abs(n-2)
w₁₂(x, y, z) = n*(n-1)*g(x, y)^abs(n-2)
w₂₂(x, y, z) = n*(n-1)*g(x, y)^abs(n-2)

w₁₁₁(x, y, z) = n*(n-1)*(n-2)*g(x, y)^abs(n-3)
w₁₁₂(x, y, z) = n*(n-1)*(n-2)*g(x, y)^abs(n-3)
w₁₂₂(x, y, z) = n*(n-1)*(n-2)*g(x, y)^abs(n-3)
w₂₂₂(x, y, z) = n*(n-1)*(n-2)*g(x, y)^abs(n-3)

w₁₁₁₁(x, y, z) = n*(n-1)*(n-2)*(n-3)*g(x, y)^abs(n-4)
w₁₁₁₂(x, y, z) = n*(n-1)*(n-2)*(n-3)*g(x, y)^abs(n-4)
w₁₁₂₂(x, y, z) = n*(n-1)*(n-2)*(n-3)*g(x, y)^abs(n-4)
w₁₂₂₂(x, y, z) = n*(n-1)*(n-2)*(n-3)*g(x, y)^abs(n-4)
w₂₂₂₂(x, y, z) = n*(n-1)*(n-2)*(n-3)*g(x, y)^abs(n-4)

κ₁₁(x,y,z) = -w₁₁(x,y,z)
κ₁₂(x, y, z) = -w₁₂(x, y, z)
κ₂₂(x, y, z) = -w₂₂(x, y, z)

M₁₁(x, y, z) = -D*(w₁₁(x, y, z) + ν*w₂₂(x, y, z))
M₂₂(x, y, z) = -D*(ν*w₁₁(x, y, z) + w₂₂(x, y, z))
M₁₂(x, y, z) = -D*(1.0 - ν)*w₁₂(x, y, z)

M₁₁₁(x, y, z) = -D*(w₁₁₁(x, y, z) + ν*w₁₂₂(x, y, z))
M₁₁₂(x, y, z) = -D*(w₁₁₂(x, y, z) + ν*w₂₂₂(x, y, z))
M₂₂₁(x, y, z) = -D*(ν*w₁₁₁(x, y, z) + w₁₂₂(x, y, z))
M₂₂₂(x, y, z) = -D*(ν*w₁₁₂(x, y, z) + w₂₂₂(x, y, z))
M₁₂₁(x, y, z) = -D*(1.0 - ν)*w₁₁₂(x, y, z)
M₁₂₂(x, y, z) = -D*(1.0 - ν)*w₁₂₂(x, y, z)

M₁₁₁₁(x, y, z) = -D*(w₁₁₁₁(x, y, z) + ν*w₁₁₂₂(x, y, z))
M₂₂₂₂(x, y, z) = -D*(ν*w₁₁₂₂(x, y, z) + w₂₂₂₂(x, y, z))
M₁₂₁₂(x, y, z) = -D*(1.0 - ν)*w₁₁₂₂(x, y, z)

q̄(x, y, z) = -(M₁₁₁₁(x, y, z)+2M₁₂₁₂(x, y, z)+M₂₂₂₂(x, y, z))
θ(x, y, z, n₁, n₂) = w₁(x, y, z)*n₁ + w₂(x, y, z)*n₂

function V̄(x, y, z, n₁, n₂)
    s₁ = -n₂; s₂ = n₁
    return (M₁₁₁(x, y, z) + M₁₂₂(x, y, z))*n₁ +
           (M₁₂₁(x, y, z) + M₂₂₂(x, y, z))*n₂ +
            M₁₁₁(x, y, z)*n₁*s₁*s₁ + M₁₁₂(x, y, z)*n₁*s₁*s₂ +
            M₂₂₁(x, y, z)*n₂*s₂*s₁ + M₂₂₂(x, y, z)*n₂*s₂*s₂ +
            M₁₂₁(x, y, z)*(n₁*s₂ + n₂*s₁)*s₁ +
            M₁₂₂(x, y, z)*(n₁*s₂ + n₂*s₁)*s₂
end

M̄(x, y, z, n₁, n₂) = M₁₁(x, y, z)*n₁*n₁ + M₂₂(x, y, z)*n₂*n₂ +
                     2.0*M₁₂(x, y, z)*n₁*n₂
ΔM̄(x,y,z) = 0.0


# 5-point Gauss-Lobatto, flat local coordinates
const ξq = [-1.0, -sqrt(3.0/7.0), 0.0, sqrt(3.0/7.0), 1.0]
const wq = [1.0/10.0, 49.0/90.0, 32.0/45.0, 49.0/90.0, 1.0/10.0]
const GL5 = (Float64[
    -1.0, 0.0, 0.0,
    -sqrt(3.0/7.0), 0.0, 0.0,
     0.0, 0.0, 0.0,
     sqrt(3.0/7.0), 0.0, 0.0,
     1.0, 0.0, 0.0,
], wq)

# ---------------------------------------------------------------------------
# Exact M monomial coefficients, in the PiecewisePolynomial basis
# ---------------------------------------------------------------------------
function exact_M_coefficients(::Val{0})
    return zeros(1), zeros(1), zeros(1)       # constant w: M = 0
end

function exact_M_coefficients(::Val{1})
    return zeros(1), zeros(1), zeros(1)       # linear w: M = 0
end

function exact_M_coefficients(::Val{2})
    c = -D * 2 * (1.0 + ν)
    d = -D * 2 * (1.0 - ν)
    return [c], [c], [d]                      # constant M
end

function exact_M_coefficients(::Val{3})
    c = -D * 3 * 2 * (1.0 + ν)
    d = -D * 3 * 2 * (1.0 - ν)
    return [c, c, c], [c, c, c], [d, d, d]   # basis 1, x, y
end

function exact_M_coefficients(::Val{4})
    c = -D * 4 * 3 * (1.0 + ν)
    d = -D * 4 * 3 * (1.0 - ν)
    shape = [1.0, 2.0, 2.0, 1.0, 2.0, 1.0]    # 1, x, y, x^2, xy, y^2
    return c .* shape, c .* shape, d .* shape
end

# ---------------------------------------------------------------------------
# Surface plot of the signed nodal residual R^w on the actual mesh nodes.
# The scattered mesh is drawn as a triangulated surface with mesh!.
# ---------------------------------------------------------------------------
function plot_rw_surface(x, y, Rw, faces; filename)
    verts = [CairoMakie.Point3f(x[i], y[i], Rw[i]) for i in 1:length(x)]
    tris = [CairoMakie.TriangleFace(Int(f[1]), Int(f[2]), Int(f[3])) for f in faces]

    # Same fixed residual range as test_gauss_consistency.jl.
    zlo, zhi = -149.0, 149.0

    fig = Figure(size = (1000, 750))
    ax = Axis3(fig[1, 1];
               xlabel = "x", ylabel = "y", zlabel = "",
               xlabelfont = "Times New Roman Italic",
               ylabelfont = "Times New Roman Italic",
               zlabelfont = "Times New Roman Italic",
               limits = ((0.0, 1.0), (0.0, 1.0), (zlo, zhi)),
               zticks = [-100.0, 0.0, 100.0],
               elevation = deg2rad(40),
               azimuth = deg2rad(-135),
               perspectiveness = 0.55,
               aspect = (1.0, 1.0, 0.6),
               xlabelsize = 26, ylabelsize = 26, zlabelsize = 26,
               xticklabelsize = 22, yticklabelsize = 22, zticklabelsize = 22)
    mesh!(ax, verts, tris; color = Rw, colormap = :balance,
          colorrange = (-150.0, 150.0), shading = true)
    Colorbar(fig[1, 2], limits = (-150.0, 150.0), colormap = :balance,
             width = 30, height = Relative(0.5),
             ticklabelsize = 22, labelsize = 26)
    mkpath(dirname(filename))
    save(filename, fig)
    println("  surface saved to $(relpath(filename, @__DIR__))")
    return fig
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
p in (2, 3, 4) || error("p must be 2, 3 or 4")
0 <= n <= p || error("the test order must satisfy 0 <= n <= p")
q >= max(0, n - 2) || error("the M basis degree must satisfy q >= n - 2")

gmsh.initialize()
gmsh.open(meshfile)
entities = getPhysicalGroups()
nodes = get𝑿ᵢ()

# Triangles for the nodal residual surface.
elemTypes, _, elemNodeTags = gmsh.model.mesh.getElements(2, -1)
faces = Tuple{Int,Int,Int}[]
for (etype, ntags) in zip(elemTypes, elemNodeTags)
    etype == 2 || continue          # 3-node triangle
    for e in 1:3:length(ntags)
        push!(faces, (Int(ntags[e]), Int(ntags[e+1]), Int(ntags[e+2])))
    end
end

nʷ = length(nodes)
area = (maximum(nodes.x) - minimum(nodes.x)) *
       (maximum(nodes.y) - minimum(nodes.y))
h = sqrt(area / nʷ)
support = (p == 2 ? 3.0 : p == 3 ? 4.0 : 5.0) * h
push!(nodes, :s₁ => fill(support, nʷ),
             :s₂ => fill(support, nʷ),
             :s₃ => fill(support, nʷ))

sp = RegularGrid(nodes.x, nodes.y, nodes.z; n = 3, γ = 5)
wtype = p == 2 ? ReproducingKernel{:Quadratic2D,:□,:QuinticSpline} :
        p == 3 ? ReproducingKernel{:Cubic2D,:□,:QuinticSpline} :
                 ReproducingKernel{:Quartic2D,:□,:QuinticSpline}
mtype = q == 0 ? PiecewisePolynomial{:Constant} :
        q == 1 ? PiecewisePolynomial{:Linear2D} :
                 PiecewisePolynomial{:Quadratic2D}

elements_w = getElements(nodes, entities["Ω"], wtype, order, sp)
prescribe!(elements_w, :q => q̄)
set∇²𝝭!(elements_w)

elements_m = getPiecewiseElements(entities["Ω"], mtype, order)
if q <= 1
    set∇𝝭!(elements_m)
    push!(elements_m, :∂²𝝭∂x², :∂²𝝭∂y², :∂²𝝭∂x∂y)
else
    set∇²𝝭!(elements_m)
end

nₑ = length(elements_m)
nM = get𝑛𝑝(elements_m[1])
nᵐ = nₑ * nM

kᵐᵐ = zeros(3*nᵐ, 3*nᵐ)
kᵐʷ = zeros(3*nᵐ, nʷ)
fʷ = zeros(nʷ)
fᵐ = zeros(3*nᵐ)

𝑎ᵐᵐ = ∫MC⁻¹MdΩ => elements_m
𝑎ᵐʷ = ∫∇²MwdΩ => (elements_m, elements_w)
𝑓ʷ  = ∫wqdΩ => elements_w

𝑎ᵐᵐ(kᵐᵐ; D = D, ν = ν)
𝑎ᵐʷ(kᵐʷ)
𝑓ʷ(fʷ)

# Γ_C blocks
elements_Γ_w = getElements(nodes, entities["Γ"], wtype, GL5, sp; normal = true)
elements_Γ_m = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], mtype, GL5)
set∇𝝭!(elements_Γ_w)
set∇𝝭!(elements_Γ_m)
𝑎ᵐʷ = [
    ∫VwdΓ  => (elements_Γ_m, elements_Γ_w),
    ∫MθdΓ => (elements_Γ_m, elements_Γ_w),
    wP    => (elements_Γ_m, elements_Γ_w),
]
𝑎ᵐʷ(kᵐʷ)

# ---------------------------------------------------------------------------
# Natural boundary data on Γᵍ: ∫wVdΓ, ∫θₙMₙₙdΓ
# ---------------------------------------------------------------------------
elements_Γ¹_w = getElements(nodes, entities["Γ¹"], wtype, GL5, sp; normal = true)
elements_Γ²_w = getElements(nodes, entities["Γ²"], wtype, GL5, sp; normal = true)
elements_Γ³_w = getElements(nodes, entities["Γ³"], wtype, GL5, sp; normal = true)
elements_Γ⁴_w = getElements(nodes, entities["Γ⁴"], wtype, GL5, sp; normal = true)
for e in (elements_Γ¹_w, elements_Γ²_w, elements_Γ³_w, elements_Γ⁴_w)
    prescribe!(e, :V => V̄, :M => M̄, :ΔM => ΔM̄)
    set∇𝝭!(e)
end
elements_Γᵍ_w = elements_Γ¹_w ∪ elements_Γ²_w ∪ elements_Γ³_w ∪ elements_Γ⁴_w
𝑓ʷ = [
    ∫wVdΓ      => elements_Γᵍ_w,
    ∫θₙMₙₙdΓ  => elements_Γᵍ_w,
]
# 𝑓ʷ(fʷ)

# ---------------------------------------------------------------------------
# Essential boundary data on Γᵍ (switch: uncomment this block and comment
# the natural block above)
# ---------------------------------------------------------------------------
elements_Γ¹_w = getElements(nodes, entities["Γ¹"], wtype, GL5, sp; normal = true)
elements_Γ²_w = getElements(nodes, entities["Γ²"], wtype, GL5, sp; normal = true)
elements_Γ³_w = getElements(nodes, entities["Γ³"], wtype, GL5, sp; normal = true)
elements_Γ⁴_w = getElements(nodes, entities["Γ⁴"], wtype, GL5, sp; normal = true)
elements_Γ¹_m = getElements(entities["Γ¹"], entities["Γ"], elements_Γ_m)
elements_Γ²_m = getElements(entities["Γ²"], entities["Γ"], elements_Γ_m)
elements_Γ³_m = getElements(entities["Γ³"], entities["Γ"], elements_Γ_m)
elements_Γ⁴_m = getElements(entities["Γ⁴"], entities["Γ"], elements_Γ_m)
for e in (elements_Γ¹_w, elements_Γ²_w, elements_Γ³_w, elements_Γ⁴_w)
    prescribe!(e, :w => w, :θ => θ)
    set∇𝝭!(e)
end
elements_Γᵍ_w = elements_Γ¹_w ∪ elements_Γ²_w ∪ elements_Γ³_w ∪ elements_Γ⁴_w
elements_Γᵍ_m = elements_Γ¹_m ∪ elements_Γ²_m ∪ elements_Γ³_m ∪ elements_Γ⁴_m
𝑎ᵐʷ = [
    ∫VwdΓ  => (elements_Γᵍ_m, elements_Γᵍ_w),
    ∫MθdΓ => (elements_Γᵍ_m, elements_Γᵍ_w),
    wP    => (elements_Γᵍ_m, elements_Γᵍ_w),
]
𝑎ᵐʷ(kᵐʷ, fᵐ)

# ---------------------------------------------------------------------------
# Exact mixed vector
# ---------------------------------------------------------------------------
dʷ = [w(nodes.x[i], nodes.y[i], 0.0) for i in 1:nʷ]

dᵐ = zeros(3*nᵐ)
M11, M22, M12 = exact_M_coefficients(Val(n))
pad_to(c, n) = length(c) == n ? c : [c; zeros(n - length(c))]
M11 = pad_to(M11, nM)
M22 = pad_to(M22, nM)
M12 = pad_to(M12, nM)
for c in 1:nₑ
    base = (c-1)*nM
    for j in 1:nM
        I = base + j
        dᵐ[3I-2] = M11[j]
        dᵐ[3I-1] = M22[j]
        dᵐ[3I]   = M12[j]
    end
end

Rᵐ = kᵐᵐ*dᵐ + kᵐʷ*dʷ - fᵐ
Rʷ = kᵐʷ'*dᵐ - fʷ

@printf("RKGSI variational-consistency residual\n")
@printf("  p = %d, n = %d, cells = %d, w DOFs = %d, M DOFs = %d\n", p, n, nₑ, nʷ, 3*nᵐ)
@printf("  ||fᵐ||             = %.6e\n", norm(fᵐ))
@printf("  ||Rᵐ||             = %.6e\n", norm(Rᵐ))
@printf("  ||Rᵐ||/max(||fᵐ||,1) = %.6e\n", norm(Rᵐ)/max(norm(fᵐ), 1.0))
@printf("  ||fʷ||             = %.6e\n", norm(fʷ))
@printf("  ||Rʷ||             = %.6e\n", norm(Rʷ))
@printf("  ||Rʷ||/max(||fʷ||,1) = %.6e\n", norm(Rʷ)/max(norm(fʷ), 1.0))
# @printf("  ||R||/max(||f||,1) = %.6e\n", norm(R)/max(norm(f), 1.0))
@printf("  Rʷ range = [%.6e, %.6e]\n", extrema(Rʷ)...)

println()
println("Plotting R^w surface ...")
plot_rw_surface(collect(nodes.x), collect(nodes.y), Rʷ, faces;
                filename = plotfile)

gmsh.finalize()
