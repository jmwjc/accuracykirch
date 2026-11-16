# ---------------------------------------------------------------------------
# patchtest_gauss.jl
#
# Standard Galerkin meshfree patch test for the Kirchhoff plate, using
# Gauss integration only (GI).  This is the baseline against which RKGSI will
# be compared.
#
# Manufactured solution on the unit square:
#
#     w(x,y) = (x + y + 1)^p
#
# with the exact normal derivative ∂w/∂n and the exact load q = D Δ²w.
# Essential boundary conditions (w and ∂w/∂n) are imposed with Nitsche's
# method.  The boundary integrals are assembled one Gmsh edge at a time,
# because ApproxOperator v0.2.0 duplicates meshfree elements when a physical
# group contains several entities.
#
# Run from the project root:
#
#     julia program/patchtest_gauss.jl
#
# The program uses the non-uniform mesh program/msh/patchtest_nonuniform.msh
# and prints the relative L2, H1 and He errors for several Gauss orders.
#
# Expected behaviour:
#   * RK shape functions are rational, and their support boundaries cut the
#     integration cells.  Gauss integration is therefore not exact even for a
#     polynomial solution.  The non-zero error is the integration-consistency
#     error of GI; this is exactly the quantity RKGSI is designed to remove.
#   * On a fixed non-uniform mesh, increasing the Gauss order reveals the
#     integration-consistency floor of GI.
# ---------------------------------------------------------------------------

using ApproxOperator
import ApproxOperator: AbstractElement
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.KirchhoffPlate: ∫κᵢⱼMᵢⱼdΩ, ∫wqdΩ, ∫VgdΓ, ∫MₙₙθdΓ, ΔMₙₛg, L₂, H₁, Hₑ
import Gmsh: gmsh

using LinearAlgebra
using Printf

# ---------------------------------------------------------------------------
# User parameters
# ---------------------------------------------------------------------------
const p           = 2                 # RK degree: 2, 3 or 4
const n           = p                 # exact solution degree: w = (x+y+1)^n
const meshfile    = joinpath(@__DIR__, "msh", "patchtest_nonuniform.msh")
const orders      = [2, 4, 6, 8]      # Gauss orders
const D           = 2.0e4             # bending stiffness
const ν           = 0.3               # Poisson ratio
const β1factor    = 1.0e3             # β1 = factor * D / h^3
const β2factor    = 1.0e3             # β2 = factor * D / h
const βcfactor    = 1.0e3             # βc = factor * D / h^2  (corner displacement)
const post_order  = 10                # quadrature order for error evaluation

# ---------------------------------------------------------------------------
# Manufactured solution
# ---------------------------------------------------------------------------
@inline f(x, y) = x + y + 1.0

w_exact(x, y, z)   = f(x, y)^n
wx_exact(x, y, z)  = n * f(x, y)^(n - 1)
wy_exact(x, y, z)  = n * f(x, y)^(n - 1)
wxx_exact(x, y, z) = n * (n - 1) * f(x, y)^(n - 2)
wxy_exact(x, y, z) = n * (n - 1) * f(x, y)^(n - 2)
wyy_exact(x, y, z) = n * (n - 1) * f(x, y)^(n - 2)

# q = D Δ²w.  For n < 4 the biharmonic is zero.
q_exact(x, y, z) = n >= 4 ?
    D * 4.0 * n * (n - 1) * (n - 2) * (n - 3) * f(x, y)^(n - 4) : 0.0

# Normal derivative on a boundary with outward normal (n1,n2).
θ_exact(x, y, z, n1, n2) = n * f(x, y)^(n - 1) * (n1 + n2)

# 5-point Gauss-Lobatto rule, used only for the corner (point) term.
const ξq = [-1.0, -sqrt(3.0/7.0), 0.0, sqrt(3.0/7.0), 1.0]
const wq = [1.0/10.0, 49.0/90.0, 32.0/45.0, 49.0/90.0, 1.0/10.0]
const GL5 = (Float64[
    -1.0, 0.0, 0.0,
    -sqrt(3.0/7.0), 0.0, 0.0,
     0.0, 0.0, 0.0,
     sqrt(3.0/7.0), 0.0, 0.0,
     1.0, 0.0, 0.0,
], wq)

# Fill the corner-jump coefficients expected by KirchhoffPlate.ΔMₙₛg on a
# meshfree boundary element.  The package's meshfree importer provides ξ.Δ
# and the normals; ΔMₙₛg additionally reads the endpoint combinations
# ξ.Δn₁s₁, ξ.Δn₁s₂n₂s₁, ξ.Δn₂s₂.
function add_corner_delta!(Γ)
    isempty(Γ) && return Γ
    data = getfield(Γ[1].𝓖[1], :data)
    data[:Δn₁s₁] = (1, Float64[])
    data[:Δn₁s₂n₂s₁] = (1, Float64[])
    data[:Δn₂s₂] = (1, Float64[])
    for ξ in Γ[1].𝓖
        Δ = ξ.Δ
        n₁ = ξ.n₁; n₂ = ξ.n₂
        s₁ = -n₂; s₂ = n₁
        push!(data[:Δn₁s₁][2], Δ*n₁*s₁)
        push!(data[:Δn₁s₂n₂s₁][2], Δ*(n₁*s₂ + n₂*s₁))
        push!(data[:Δn₂s₂][2], Δ*n₂*s₂)
    end
    return Γ
end

# ---------------------------------------------------------------------------
# One patch-test run
# ---------------------------------------------------------------------------
function run_patchtest(order::Int)
    p in (2, 3, 4) || error("p must be 2, 3 or 4.")
    0 <= n <= p || error("the exact solution degree must satisfy 0 <= n <= p.")

    # Suppress the Gmsh and importer summaries; only the result table is shown.
    return redirect_stdout(devnull) do
        gmsh.initialize()
        gmsh.open(meshfile)
        entities = getPhysicalGroups()
        nodes = get𝑿ᵢ()

        x, y, z = nodes.x, nodes.y, nodes.z
        nₚ = length(nodes)
        area = (maximum(x) - minimum(x)) * (maximum(y) - minimum(y))
        h = sqrt(area / nₚ)

        # The graded mesh has coarser regions, where the 3.5 h cubic support
        # can make the moment matrix nearly singular.  Use a slightly larger
        # normalized support than on the uniform mesh.
        support = (p == 2 ? 3.0 : p == 3 ? 4.0 : 5.0) * h
        push!(nodes, :s₁ => fill(support, nₚ),
                     :s₂ => fill(support, nₚ),
                     :s₃ => fill(support, nₚ))

        sp = RegularGrid(x, y, z; n = 3, γ = 5)
        type = p == 2 ? ReproducingKernel{:Quadratic2D,:□,:QuinticSpline} :
               p == 3 ? ReproducingKernel{:Cubic2D,:□,:QuinticSpline} :
                        ReproducingKernel{:Quartic2D,:□,:QuinticSpline}

        # --- domain: standard Gauss integration ----------------------------
        Ω = getElements(nodes, entities["Ω"], type, order, sp)
        prescribe!(Ω, :q => q_exact)
        set∇²𝝭!(Ω)

        k = zeros(nₚ, nₚ)
        f = zeros(nₚ)

        (∫κᵢⱼMᵢⱼdΩ => Ω)(k; D = D, ν = ν)
        (∫wqdΩ       => Ω)(f)

        ngΩ = length(Ω[1].𝓖)      # Gauss points per domain cell

        # --- boundary: Nitsche, one Gmsh edge at a time --------------------
        β1 = β1factor * D / h^3
        β2 = β2factor * D / h
        βc = βcfactor * D / h^2
        ngΓ = 0

        for tag in entities["Γᵍ"].second
            Γ = getElements(nodes, 1 => [tag], type, order, sp; normal = true)
            ngΓ = length(Γ[1].𝓖)  # Gauss points per boundary cell
            prescribe!(Γ, :g => w_exact, :θ => θ_exact)
            set∇̂³𝝭!(Γ)

            # Package Nitsche operators on the edge:
            #   displacement : ∫VgdΓ       (needs set∇̂³𝝭!)
            #   rotation     : ∫MₙₙθdΓ     (needs set∇²𝝭!, included in ∇̂³)
            (∫VgdΓ   => Γ)(k, f; D = D, ν = ν, α = β1)
            (∫MₙₙθdΓ => Γ)(k, f; D = D, ν = ν, α = β2)

            # Corner displacement term: ΔMₙₛg is a point contribution, so
            # use the Lobatto rule, whose endpoints carry ξ.Δ ≠ 0.
            Γc = getElements(nodes, 1 => [tag], type, GL5, sp;
                             normal = true, delta = true)
            add_corner_delta!(Γc)
            prescribe!(Γc, :g => w_exact)
            set∇²𝝭!(Γc)
            (ΔMₙₛg => Γc)(k, f; D = D, ν = ν, α = βc)
        end

        # Residual of the discrete equations for the exact nodal vector.
        # This is the direct integration-consistency measure of the patch test.
        d_exact = [w_exact(nodes.x[i], nodes.y[i], 0.0) for i in 1:nₚ]
        resid = norm(k * d_exact - f) / norm(f)

        d = k \ f

        # --- error evaluation with the package operators --------------------
        # L₂, H₁ and Hₑ read the solution from xᵢ.d, so push d onto the nodes.
        push!(nodes, :d => d)
        Ωp = getElements(nodes, entities["Ω"], type, post_order, sp)
        set∇²𝝭!(Ωp)
        prescribe!(Ωp,
                   :w        => w_exact,
                   :∂w∂x     => wx_exact,
                   :∂w∂y     => wy_exact,
                   :∂²w∂x²   => wxx_exact,
                   :∂²w∂x∂y  => wxy_exact,
                   :∂²w∂y²   => wyy_exact)

        L2 = L₂(Ωp)
        H1, _ = H₁(Ωp)
        He = Hₑ(Ωp; D = D, ν = ν)

        gmsh.finalize()
        (ngΩ, ngΓ, L2, H1, He, resid)
    end
end

# ---------------------------------------------------------------------------
# Main: Gauss-order sweep on the fixed non-uniform mesh
# ---------------------------------------------------------------------------
println("Gauss-integration Galerkin meshfree patch test (non-uniform mesh)")
@printf("  p = %d, n = %d, mesh = %s\n", p, n, basename(meshfile))
@printf("  D = %.3e, nu = %.2f\n", D, ν)
@printf("  Nitsche: beta1 = %.1e * D/h^3, beta2 = %.1e * D/h, betac = %.1e * D/h^2\n",
        β1factor, β2factor, βcfactor)
@printf("  support = %.1f h, post-processing order = %d\n",
        p == 2 ? 3.0 : p == 3 ? 4.0 : 5.0, post_order)
println()
println("  n_gauss(Ω)  n_gauss(Γ)        L2 error        H1 error        He error       residual")
println("  -----------------------------------------------------------------------------------------")

for order in orders
    ngΩ, ngΓ, L2, H1, He, resid = run_patchtest(order)
    @printf("  %10d  %10d   %12.5e   %12.5e   %12.5e   %12.5e\n",
            ngΩ, ngΓ, L2, H1, He, resid)
end

println()
println("Notes:")
println("  * Errors are computed with the package operators")
println("    KirchhoffPlate.L₂, H₁ and Hₑ.  Hₑ is the relative energy norm")
println("    (D-weighted curvature error), the package's H²-equivalent measure.")
println("  * GI does not reach round-off because RK shape functions are rational")
println("    and their support boundaries cut the integration cells; the residual")
println("    is the integration-consistency error that RKGSI targets.")
println("  * Columns n_gauss(Ω) and n_gauss(Γ) are the Gauss-point counts")
println("    per triangular cell and per boundary edge, respectively.")
println("  * p = 2, 3 and 4 use set∇̂³𝝭! for the boundary Kirchhoff-shear term;")
println("    on this graded mesh the normalized supports are 3.0 h, 4.0 h and 5.0 h.")
