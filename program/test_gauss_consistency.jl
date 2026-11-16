# ---------------------------------------------------------------------------
# test_gauss_consistency.jl
#
# Variational-consistency test for the standard Gauss-integrated Galerkin
# discretization of the Kirchhoff plate.
#
# The exact polynomial
#
#     w = (x + y + 1)^n
#
# is inserted into the discrete equations.  For a variationally consistent
# discretization the residual
#
#     R = K d_exact - f
#
# must vanish.  With Gauss integration of the RK shape functions it does not:
# the residual is the integration-consistency error of GI, which is the
# quantity RKGSI is designed to remove.
#
# The Nitsche boundary terms and the corner displacement term are the same as
# in patchtest_gauss.jl:
#
#     ∫VgdΓ,  ∫MₙₙθdΓ,  ΔMₙₛg
#
# Run from the project root:
#
#     julia program/test_gauss_consistency.jl
# ---------------------------------------------------------------------------

using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.KirchhoffPlate:
    ∫κᵢⱼMᵢⱼdΩ, ∫wqdΩ, ∫VgdΓ, ∫MₙₙθdΓ, ∫∂ₙwθdΓ, ΔMₙₛg
import ApproxOperator.MindlinPlate: ∫αwwdΓ
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
const p           = 3                 # RK degree: 2, 3 or 4
const n           = 3                 # exact solution degree: w = (x+y+1)^n
const meshfile    = joinpath(@__DIR__, "msh", "patchtest_nonuniform.msh")
const orders      = [2, 4, 6, 8]      # Gauss orders in the domain
const D           = 1.0
const ν           = 0.3
const boundary_method = :nitsche   # :nitsche or :penalty
const β1factor    = 1.0e3             # β1 = factor * D / h^3
const β2factor    = 1.0e3             # β2 = factor * D / h
const βcfactor    = 1.0e3             # βc = factor * D / h^2  (corner w)
const plot_order  = 8                 # Gauss order used for the surface plot
const plotfile    = joinpath(@__DIR__, "figures", "gauss_consistency_surface.png")

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

q_exact(x, y, z) = n >= 4 ?
    D * 4.0 * n * (n - 1) * (n - 2) * (n - 3) * f(x, y)^(n - 4) : 0.0

θ_exact(x, y, z, n₁, n₂) = n * f(x, y)^(n - 1) * (n₁ + n₂)

# ---------------------------------------------------------------------------
# 5-point Gauss-Lobatto rule, used for the corner (point) term only
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
# Fill the corner-jump coefficients expected by KirchhoffPlate.ΔMₙₛg.
# The meshfree importer provides ξ.Δ and the edge normals; ΔMₙₛg additionally
# reads the endpoint combinations ξ.Δn₁s₁, ξ.Δn₁s₂n₂s₁, ξ.Δn₂s₂.
# ---------------------------------------------------------------------------
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
# Assemble the Gauss-Galerkin system for one Gauss order and evaluate the
# residual of the exact nodal vector.
# ---------------------------------------------------------------------------
function gauss_residual(order::Int)
    p in (2, 3, 4) || error("p must be 2, 3 or 4.")
    0 <= n <= p || error("the exact solution degree must satisfy 0 <= n <= p.")
    boundary_method in (:nitsche, :penalty) ||
        error("boundary_method must be :nitsche or :penalty.")

    return redirect_stdout(devnull) do
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

        x, y, z = nodes.x, nodes.y, nodes.z
        nₚ = length(nodes)
        area = (maximum(x) - minimum(x)) * (maximum(y) - minimum(y))
        h = sqrt(area / nₚ)

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

        ngΩ = length(Ω[1].𝓖)

        # --- boundary: Nitsche edge terms + corner displacement -----------
        β1 = β1factor * D / h^3
        β2 = β2factor * D / h
        βc = βcfactor * D / h^2
        ngΓ = 0

        for tag in entities["Γᵍ"].second
            # Edge Nitsche terms, integrated with the requested Gauss order.
            Γ = getElements(nodes, 1 => [tag], type, order, sp; normal = true)
            ngΓ = length(Γ[1].𝓖)
            prescribe!(Γ, :g => w_exact, :θ => θ_exact)
            set∇̂³𝝭!(Γ)

            if boundary_method == :nitsche
                # Consistency + penalty (Nitsche).
                (∫VgdΓ   => Γ)(k, f; D = D, ν = ν, α = β1)
                (∫MₙₙθdΓ => Γ)(k, f; D = D, ν = ν, α = β2)
            else
                # Penalty only.
                (∫αwwdΓ   => Γ)(k, f; α = β1)
                (∫∂ₙwθdΓ => Γ)(k, f; α = β2)
            end

            # Corner displacement term: ΔMₙₛg is a point contribution, so
            # use the Lobatto rule, whose endpoints carry ξ.Δ ≠ 0.
            Γc = getElements(nodes, 1 => [tag], type, GL5, sp;
                             normal = true, delta = true)
            add_corner_delta!(Γc)
            prescribe!(Γc, :g => w_exact)
            set∇²𝝭!(Γc)
            if boundary_method == :nitsche
                (ΔMₙₛg => Γc)(k, f; D = D, ν = ν, α = βc)
            else
                # D = 0 removes the consistency part, leaving α_c w δw.
                (ΔMₙₛg => Γc)(k, f; D = 0.0, ν = ν, α = βc)
            end
        end

        # Exact nodal vector and residual.
        d_exact = [w_exact(nodes.x[i], nodes.y[i], 0.0) for i in 1:nₚ]
        R = k * d_exact - f

        gmsh.finalize()
        (ngΩ, ngΓ, norm(f), norm(R), collect(x), collect(y), R, faces)
    end
end

# ---------------------------------------------------------------------------
# Surface plot of the signed nodal residual R_I on the actual mesh nodes.
#
# The nodes are scattered, so Makie's regular-grid `surface` recipe cannot be
# used directly.  The triangulated mesh from Gmsh is therefore drawn as a
# surface with `mesh!`: the vertices are exactly (x_I, y_I, R_I).
# ---------------------------------------------------------------------------
function plot_residual_surface(x, y, R, faces; order, filename)
    verts = [CairoMakie.Point3f(x[i], y[i], R[i]) for i in 1:length(x)]
    tris = [CairoMakie.TriangleFace(Int(f[1]), Int(f[2]), Int(f[3])) for f in faces]

    fig = Figure(size = (1000, 750))
    ax = Axis3(fig[1, 1];
               xlabel = "x", ylabel = "y", zlabel = "",
               xlabelfont = "Times New Roman Italic",
               ylabelfont = "Times New Roman Italic",
               zlabelfont = "Times New Roman Italic",
               limits = ((0.0, 1.0), (0.0, 1.0), (-149, 149)),
               zticks = [-100.0, 0.0, 100.0],
            #    elevation = deg2rad(32),
            #    azimuth = deg2rad(-55),
               elevation = deg2rad(40),
               azimuth = deg2rad(-135),
               perspectiveness = 0.55,
               aspect = (1.0, 1.0, 0.6),
               xlabelsize = 26, ylabelsize = 26, zlabelsize = 26,
               xticklabelsize = 22, yticklabelsize = 22, zticklabelsize = 22)
    mesh!(ax, verts, tris; color = R, colormap = :balance,
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
# Main: Gauss-order sweep
# ---------------------------------------------------------------------------
println("Gauss-integration Galerkin variational-consistency residual")
@printf("  p = %d, n = %d, mesh = %s\n", p, n, basename(meshfile))
@printf("  D = %.3e, nu = %.2f\n", D, ν)
@printf("  boundary enforcement = %s\n", boundary_method)
@printf("  Nitsche: beta1 = %.1e * D/h^3, beta2 = %.1e * D/h, betac = %.1e * D/h^2\n",
        β1factor, β2factor, βcfactor)
println()
println("  n_gauss(Ω)  n_gauss(Γ)        ||f||           ||R||        ||R||/||f||")
println("  ---------------------------------------------------------------------------------")

plot_data = nothing
for order in orders
    ngΩ, ngΓ, nf, nR, xnodes, ynodes, Rnodes, faces = gauss_residual(order)
    @printf("  %10d  %10d   %12.5e   %12.5e   %12.5e\n",
            ngΩ, ngΓ, nf, nR, nR/max(nf, 1.0))
    if order == plot_order
        global plot_data = (xnodes, ynodes, Rnodes, faces, ngΩ, ngΓ)
    end
end

if plot_data !== nothing
    println()
    println("Plotting R_I surface: method = $boundary_method, " *
            "Gauss order = $plot_order, n_gauss(Ω) = $(plot_data[5]), " *
            "n_gauss(Γ) = $(plot_data[6]), corner = 5 (Lobatto)")
    plot_residual_surface(plot_data[1], plot_data[2], plot_data[3], plot_data[4];
                          order = plot_order, filename = plotfile)
end

println()
println("Notes:")
println("  * R = K d_exact - f is the residual of the exact nodal vector.")
println("  * Gauss integration is not exact for the rational RK shape functions,")
println("    so R does not reach round-off; this is the GI consistency error.")
println("  * The corner term uses ΔMₙₛg with a 5-point Lobatto rule, whose")
println("    endpoints carry the corner sign ξ.Δ.")
