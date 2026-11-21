# ---------------------------------------------------------------------------
# patchtest_rkgsi.jl
#
# Mixed Hellinger-Reissner RKGSI patch test for the Kirchhoff plate.
#
# Unknowns
#   w : RK scalar, one DOF per mesh node
#   M : independent piecewise polynomial per integration cell,
#       three DOFs per M node (M11, M22, M12)
#
# The full saddle-point system is assembled and solved directly; no Schur
# complement is formed.
#
# Boundary operators are overloaded:
#   (aᵐ, aʷ, k)       -> cell boundary Γ_C
#   (aᵐ, aʷ, k, f)    -> essential boundary Γ_w / Γ_θ / c_w
#
# Run from the project root:
#
#     julia program/patchtest_rkgsi.jl
# ---------------------------------------------------------------------------

# macOS: Homebrew Julia, Gmsh and UMFPACK may each bundle libomp.
# Allow the duplicate runtime; see OMP Error #15.
ENV["KMP_DUPLICATE_LIB_OK"] = "TRUE"

using ApproxOperator
import ApproxOperator: AbstractElement, get𝑛𝑝
import ApproxOperator.GmshImport:
    getPhysicalGroups, get𝑿ᵢ, getElements,
    getPiecewiseElements, getPiecewiseBoundaryElements
import ApproxOperator.KirchhoffPlate:
    ∫MC⁻¹MdΩ, ∫∇²MwdΩ,
    ∫VwdΓ, ∫MθdΓ, wP,
    ∫wVdΓ, ∫θₙMₙₙdΓ, wΔMₙₛ,
    ∫wqdΩ, L₂, H₁, Hₑ
import ApproxOperator.MindlinPlate: ∫αwwdΓ
import Gmsh: gmsh

using LinearAlgebra
using Printf

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
const p          = 2                    # RK basis degree: 2, 3 or 4
const n          = 1                    # exact solution order: w = (x+y+1)^n
const q          = 0                    # M basis degree: 0 constant, 1 linear, 2 quadratic
const meshfile   = joinpath(@__DIR__, "msh", "patchtest_rkgsi.msh")
const D          = 1.0
const ν          = 0.3
const order      = p == 2 ? 4 : p == 3 ? 6 : 8
const post_order = 10

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

κ₁₁(x, y, z) = -w₁₁(x, y, z)
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

q̄(x, y, z) = -(M₁₁₁₁(x, y, z) + 2M₁₂₁₂(x, y, z) + M₂₂₂₂(x, y, z))
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
ΔM̄(x, y, z) = 0.0

# ---------------------------------------------------------------------------
# 5-point Gauss-Lobatto, flat local coordinates
# ---------------------------------------------------------------------------
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
# Main
# ---------------------------------------------------------------------------
p in (2, 3, 4) || error("p must be 2, 3 or 4")
q in (0, 1, 2) || error("q must be 0, 1 or 2")
0 <= n <= p || error("the exact solution order must satisfy 0 <= n <= p")
q >= max(0, n - 2) || error("the M basis degree must satisfy q >= n - 2")

gmsh.initialize()
gmsh.open(meshfile)
entities = getPhysicalGroups()
nodes = get𝑿ᵢ()

x, y, z = nodes.x, nodes.y, nodes.z
nʷ = length(nodes)
area = (maximum(x) - minimum(x)) * (maximum(y) - minimum(y))
h = sqrt(area / nʷ)

support = (p == 2 ? 3.0 : p == 3 ? 4.0 : 5.0) * h
αw = 1.0e3 * D / h^3
push!(nodes, :s₁ => fill(support, nʷ),
             :s₂ => fill(support, nʷ),
             :s₃ => fill(support, nʷ))

sp = RegularGrid(x, y, z; n = 3, γ = 5)
wtype = p == 2 ? ReproducingKernel{:Quadratic2D,:□,:QuinticSpline} :
        p == 3 ? ReproducingKernel{:Cubic2D,:□,:QuinticSpline} :
                 ReproducingKernel{:Quartic2D,:□,:QuinticSpline}
mtype = q == 0 ? PiecewisePolynomial{:Constant} :
        q == 1 ? PiecewisePolynomial{:Linear2D} :
                 PiecewisePolynomial{:Quadratic2D}

# domain cells from Gmsh, one per RK element
elements_w = getElements(nodes, entities["Ω"], wtype, order, sp)
prescribe!(elements_w, :q => q̄)
set∇²𝝭!(elements_w)

# independent piecewise M field on the domain
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
# Essential boundary data on Γᵍ
# (switch: uncomment this block and comment the natural block above)
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
    prescribe!(e, :w => w, :θ => θ, :g => w)
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
# Penalty stabilization of the w block for q = 0
# ---------------------------------------------------------------------------
kʷʷ = zeros(nʷ, nʷ)
if q == 0
    (∫αwwdΓ => elements_Γᵍ_w)(kʷʷ, fʷ; α = αw)
end

# ---------------------------------------------------------------------------
# Solve the full saddle-point system
# ---------------------------------------------------------------------------
d = [kᵐᵐ kᵐʷ; kᵐʷ' kʷʷ] \ [fᵐ; fʷ]
w_sol = d[3*nᵐ+1:3*nᵐ+nʷ]

# ---------------------------------------------------------------------------
# Error norms with package operators
# ---------------------------------------------------------------------------
push!(nodes, :d => w_sol)
Ωp = getElements(nodes, entities["Ω"], wtype, post_order, sp)
set∇²𝝭!(Ωp)
prescribe!(Ωp,
           :w       => w,
           :∂w∂x    => w₁,
           :∂w∂y    => w₂,
           :∂²w∂x²  => w₁₁,
           :∂²w∂x∂y => w₁₂,
           :∂²w∂y²  => w₂₂)
L2 = L₂(Ωp)
H1, _ = H₁(Ωp)
He = Hₑ(Ωp; D = D, ν = ν)

@printf("Mixed HR-RKGSI patch test\n")
@printf("  p = %d, n = %d, q = %d, cells = %d, w DOFs = %d, M DOFs = %d\n",
        p, n, q, nₑ, nʷ, 3*nᵐ)
@printf("  support = %.1f h, domain order = %d, post order = %d\n",
        p == 2 ? 3.0 : p == 3 ? 4.0 : 5.0, order, post_order)
@printf("  q = %d, penalty αw = %.3e\n", q, q == 0 ? αw : 0.0)
@printf("  L2 = %.6e   H1 = %.6e   He = %.6e\n", L2, H1, He)

gmsh.finalize()
