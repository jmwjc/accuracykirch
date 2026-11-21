# ---------------------------------------------------------------------------
# sinusoidal_gauss.jl
#
# Simply supported rectangular plate under a sinusoidal body force,
# solved with the standard Galerkin meshfree formulation using Gauss
# integration and Nitsche's method for the essential condition w = 0.
#
# This is the Gauss-integration counterpart of sinusoidal_rkgsi.jl:
#   * same geometry, material, exact solution and mesh;
#   * same RK approximation and support sizes;
#   * no independent moment field and no mixed system;
#   * w = 0 on Γ is enforced by Nitsche's method, while Mₙ = 0 is natural.
#
# Usage:
#
#     julia program/sinusoidal_gauss.jl [p] [ndiv]
#
#     p     RK basis degree: 2, 3 or 4          (default 3)
#     ndiv  number of subdivisions per edge     (default 10)
#
# Sweep mode:
#
#     julia program/sinusoidal_gauss.jl sweep [p] [outbase]
#
#     p      2, 3, 4 or all (default all)
#     ndiv   fixed by SWEEP_NDIVS = (8, 16, 32, 64)
#
# The default CSV prefix contains a timestamp, and write_results refuses to
# overwrite any existing file.  Pass an explicit new prefix if needed.
#
# Run from the project root.
# ---------------------------------------------------------------------------

# macOS: Homebrew Julia, Gmsh and LAPACK may each bundle libomp.
ENV["KMP_DUPLICATE_LIB_OK"] = "TRUE"

using ApproxOperator
import ApproxOperator: AbstractElement
import ApproxOperator.GmshImport:
    getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.KirchhoffPlate:
    ∫κᵢⱼMᵢⱼdΩ, ∫wqdΩ,
    ∫VgdΓ, ΔMₙₛg,
    L₂, H₁, Hₑ
import BenchmarkExample
import Gmsh: gmsh

using LinearAlgebra
using Printf
using Dates

# ---------------------------------------------------------------------------
# Physical and geometric parameters
# ---------------------------------------------------------------------------
const D          = 1.0
const ν          = 0.3
const L          = 1.0
const q̄₀         = 4.0*π^4
const post_order = 12

const β1factor = 1.0e3      # β1 = factor * D / h^3
const βcfactor = 1.0e3      # βc = factor * D / h^2

const α = π/L
const β = π/L
const A = q̄₀/(D*(α^2 + β^2)^2)

# ---------------------------------------------------------------------------
# Exact fields
# ---------------------------------------------------------------------------
w(x, y)      = A*sin(α*x)*sin(β*y)
w₁(x, y)     = A*α*cos(α*x)*sin(β*y)
w₂(x, y)     = A*β*sin(α*x)*cos(β*y)
w₁₁(x, y)    = -A*α^2*sin(α*x)*sin(β*y)
w₁₂(x, y)    = A*α*β*cos(α*x)*cos(β*y)
w₂₂(x, y)    = -A*β^2*sin(α*x)*sin(β*y)

q̄(x, y)      = q̄₀*sin(α*x)*sin(β*y)

# 5-point Gauss--Lobatto rule for the corner point term.
const wq = [1.0/10.0, 49.0/90.0, 32.0/45.0, 49.0/90.0, 1.0/10.0]
const GL5 = (Float64[
    -1.0, 0.0, 0.0,
    -sqrt(3.0/7.0), 0.0, 0.0,
     0.0, 0.0, 0.0,
     sqrt(3.0/7.0), 0.0, 0.0,
     1.0, 0.0, 0.0,
], wq)

# ---------------------------------------------------------------------------
# Mesh helpers (same mesh as sinusoidal_rkgsi.jl)
# ---------------------------------------------------------------------------
const mshdir = joinpath(@__DIR__, "msh")

function meshfile(ndiv::Int)
    tagL = replace(string(L), "." => "p")
    return joinpath(mshdir, "sinusoidal_L$(tagL)_$(ndiv).msh")
end

function ensure_mesh(ndiv::Int)
    mkpath(mshdir)
    file = meshfile(ndiv)
    if !isfile(file)
        generate_mesh(file, ndiv)
    end
    return file
end

function generate_mesh(file::String, ndiv::Int)
    gmsh.initialize()
    gmsh.model.add("Sinusoidal plate")

    lc = L/ndiv
    gmsh.model.geo.addPoint(0.0, 0.0, 0.0, lc, 1)
    gmsh.model.geo.addPoint(  L, 0.0, 0.0, lc, 2)
    gmsh.model.geo.addPoint(  L,   L, 0.0, lc, 3)
    gmsh.model.geo.addPoint(0.0,   L, 0.0, lc, 4)
    Γ¹ = gmsh.model.geo.addLine(1, 2, 1)
    Γ² = gmsh.model.geo.addLine(2, 3, 2)
    Γ³ = gmsh.model.geo.addLine(3, 4, 3)
    Γ⁴ = gmsh.model.geo.addLine(4, 1, 4)
    gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)
    Ω = gmsh.model.geo.addPlaneSurface([1], 1)
    gmsh.model.geo.synchronize()

    gmsh.model.addPhysicalGroup(1, [Γ¹], -1, "Γ¹")
    gmsh.model.addPhysicalGroup(1, [Γ²], -1, "Γ²")
    gmsh.model.addPhysicalGroup(1, [Γ³], -1, "Γ³")
    gmsh.model.addPhysicalGroup(1, [Γ⁴], -1, "Γ⁴")
    gmsh.model.addPhysicalGroup(2, [Ω], -1, "Ω")

    gmsh.model.mesh.setTransfiniteCurve(Γ¹, ndiv + 1)
    gmsh.model.mesh.setTransfiniteCurve(Γ², ndiv + 1)
    gmsh.model.mesh.setTransfiniteCurve(Γ³, ndiv + 1)
    gmsh.model.mesh.setTransfiniteCurve(Γ⁴, ndiv + 1)
    gmsh.model.mesh.setTransfiniteSurface(Ω)
    gmsh.model.mesh.setAlgorithm(2, Ω, 1)
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(1)

    tag = BenchmarkExample.addEdgeElements((2, 1), 1)
    gmsh.model.geo.addPhysicalGroup(1, [tag], -1, "Γ")
    gmsh.model.geo.synchronize()
    gmsh.write(file)
    gmsh.finalize()
end

# ---------------------------------------------------------------------------
# Approximation and boundary data
# ---------------------------------------------------------------------------
function wtype(p::Int)
    return p == 2 ? ReproducingKernel{:Quadratic2D,:□,:QuinticSpline} :
           p == 3 ? ReproducingKernel{:Cubic2D,:□,:QuinticSpline} :
                    ReproducingKernel{:Quartic2D,:□,:QuinticSpline}
end

function support_ratio(p::Int)
    return p == 2 ? 2.5 : p == 3 ? 3.5 : 4.5
end

# Fill the corner-jump coefficients expected by KirchhoffPlate.ΔMₙₛg.
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
# One run
# ---------------------------------------------------------------------------
function run_case(p::Int, ndiv::Int)
    p in (2, 3, 4) || error("p must be 2, 3 or 4")
    ndiv > 1       || error("ndiv must be larger than one")

    file = meshfile(ndiv)
    isfile(file) || error("mesh not found: $file; call ensure_mesh($ndiv) first")

    gmsh.initialize()
    gmsh.open(file)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()

    x, y, z = nodes.x, nodes.y, nodes.z
    nʷ = length(nodes)
    h = 1.0/ndiv
    support = support_ratio(p)*h
    push!(nodes, :s₁ => fill(support, nʷ),
                 :s₂ => fill(support, nʷ),
                 :s₃ => fill(support, nʷ))

    sp = RegularGrid(x, y, z; n = 3, γ = 5)
    # order = p == 2 ? 4 : p == 3 ? 6 : 8
    order = p == 2 ? 2 : p == 3 ? 3 : 4

    # --- domain: Gauss integration ---------------------------------------
    elements = getElements(nodes, entities["Ω"], wtype(p), order, sp)
    prescribe!(elements, :q => q̄)
    set∇²𝝭!(elements)

    k = zeros(nʷ, nʷ)
    f = zeros(nʷ)

    (∫κᵢⱼMᵢⱼdΩ => elements)(k; D = D, ν = ν)
    (∫wqdΩ       => elements)(f)

    # --- simply supported boundary: Nitsche for w = 0 --------------------
    β1 = β1factor * D / h^3
    βc = βcfactor * D / h^2

    for key in ("Γ¹", "Γ²", "Γ³", "Γ⁴")
        Γ = getElements(nodes, entities[key], wtype(p), order, sp; normal = true)
        prescribe!(Γ, :g => w)
        set∇̂³𝝭!(Γ)
        (∫VgdΓ => Γ)(k, f; D = D, ν = ν, α = β1)

        Γc = getElements(nodes, entities[key], wtype(p), GL5, sp;
                         normal = true, delta = true)
        add_corner_delta!(Γc)
        prescribe!(Γc, :g => w)
        set∇²𝝭!(Γc)
        (ΔMₙₛg => Γc)(k, f; D = D, ν = ν, α = βc)
    end

    d = k \ f

    # --- errors -----------------------------------------------------------
    push!(nodes, :d => d)
    Ωp = getElements(nodes, entities["Ω"], wtype(p), post_order, sp)
    set∇²𝝭!(Ωp)
    prescribe!(Ωp,
               :w        => w,
               :∂w∂x     => w₁,
               :∂w∂y     => w₂,
               :∂²w∂x²   => w₁₁,
               :∂²w∂x∂y  => w₁₂,
               :∂²w∂y²   => w₂₂)

    L2 = L₂(Ωp)
    H1, _ = H₁(Ωp)
    He = Hₑ(Ωp; D = D, ν = ν)

    gmsh.finalize()

    return (; p, ndiv, h, nʷ, L2, H1, H2 = He, He)
end

# ---------------------------------------------------------------------------
# CSV output (no overwrite)
# ---------------------------------------------------------------------------
const SWEEP_PS    = (2, 3, 4)
const SWEEP_NDIVS = (8, 16, 32, 64)

function timestamp()
    return Dates.format(now(), "yyyymmdd_HHMMSS")
end

function default_outbase()
    return joinpath(@__DIR__, "results", "sinusoidal_gauss_" * timestamp())
end

function ensure_new_output(outbase)
    for suffix in ("", "_rates", "_config")
        file = outbase * suffix * ".csv"
        if isfile(file)
            error("refusing to overwrite existing result: $file; choose a new output prefix")
        end
    end
end

function csv_cell(x)
    s = ismissing(x) ? "" : string(x)
    if occursin(',', s) || occursin('"', s) || occursin('\n', s)
        return "\"" * replace(s, "\"" => "\"\"") * "\""
    end
    return s
end

function write_csv(file, header, cols)
    mkpath(dirname(file))
    open(file, "w") do io
        println(io, join(header, ","))
        for i in eachindex(cols[1])
            println(io, join([csv_cell(cols[j][i]) for j in eachindex(cols)], ","))
        end
    end
    return file
end

function write_results(rows, outbase = default_outbase())
    ensure_new_output(outbase)

    raw = (
        p     = Int[r.p for r in rows],
        ndiv  = Int[r.ndiv for r in rows],
        nodes = Int[r.nʷ for r in rows],
        h     = Any[r.h for r in rows],
        L2    = Any[r.L2 for r in rows],
        H1    = Any[r.H1 for r in rows],
        H2    = Any[r.H2 for r in rows],
        He    = Any[r.He for r in rows],
        valid = String[(r.L2 isa Real && r.H1 isa Real && r.He isa Real &&
                        isfinite(r.L2) && isfinite(r.H1) && isfinite(r.He)) ? "yes" : "no"
                       for r in rows],
    )

    rate_p = Int[]
    rate_norm = String[]
    r_8_16 = Any[]
    r_16_32 = Any[]
    r_32_64 = Any[]
    r_all = Any[]

    for p in unique(r.p for r in rows)
        sub = sort([r for r in rows if r.p == p], by = r -> r.ndiv)
        length(sub) == 4 || continue
        ok = all(r -> r.L2 isa Real && r.H1 isa Real && r.He isa Real &&
                       isfinite(r.L2) && isfinite(r.H1) && isfinite(r.He) &&
                       r.L2 > 0 && r.H1 > 0 && r.He > 0, sub)
        if !ok
            for name in ("L2", "H1", "H2")
                push!(rate_p, p)
                push!(rate_norm, name)
                push!(r_8_16, ""); push!(r_16_32, ""); push!(r_32_64, ""); push!(r_all, "")
            end
            continue
        end
        for (name, key) in (("L2", :L2), ("H1", :H1), ("H2", :H2))
            vals = [getproperty(r, key) for r in sub]
            push!(rate_p, p)
            push!(rate_norm, name)
            push!(r_8_16,  log(vals[1]/vals[2])/log(2))
            push!(r_16_32, log(vals[2]/vals[3])/log(2))
            push!(r_32_64, log(vals[3]/vals[4])/log(2))
            push!(r_all,   log(vals[1]/vals[4])/log(8))
        end
    end

    rates = (; p = rate_p, norm = rate_norm,
               rate_8_16 = r_8_16, rate_16_32 = r_16_32,
               rate_32_64 = r_32_64, rate_all = r_all)

    config = (
        item = ["method", "D", "nu", "L", "q0", "error norm", "support", "Nitsche"],
        value = ["Gauss integration + Nitsche",
                 string(D), string(ν), string(L), string(q̄₀),
                 "H2 = He = relative energy norm of w - w^h",
                 "p=2: 2.5h; p=3: 3.5h; p=4: 4.5h",
                 "β1 = 1e3 D/h^3, βc = 1e3 D/h^2; w=0 only"],
    )

    raw_file = write_csv(outbase * ".csv", String.(keys(raw)),
                         [Any[x for x in getproperty(raw, k)] for k in keys(raw)])
    rate_file = write_csv(outbase * "_rates.csv", String.(keys(rates)),
                          [Any[x for x in getproperty(rates, k)] for k in keys(rates)])
    config_file = write_csv(outbase * "_config.csv", String.(keys(config)),
                            [Any[x for x in getproperty(config, k)] for k in keys(config)])

    @printf("written: %s\n", raw_file)
    @printf("written: %s\n", rate_file)
    @printf("written: %s\n", config_file)
    return raw_file
end

function run_sweep(ps = SWEEP_PS, outbase = default_outbase())
    ensure_mesh(8)
    ensure_mesh(16)
    ensure_mesh(32)
    ensure_mesh(64)

    rows = NamedTuple[]
    for p in ps, ndiv in SWEEP_NDIVS
        @printf("running p = %d, ndiv = %d ...\n", p, ndiv)
        flush(stdout)
        r = try
            run_case(p, ndiv)
        catch err
            @warn "run failed" p ndiv exception = err
            nothing
        end
        if r === nothing
            push!(rows, (; p, ndiv, h = "", nʷ = 0,
                            L2 = "", H1 = "", H2 = "", He = ""))
        else
            push!(rows, r)
        end
    end

    write_results(rows, outbase)
    return rows
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
function main()
    if get(ARGS, 1, "") == "mesh"
        ensure_mesh(8)
        ensure_mesh(16)
        ensure_mesh(32)
        ensure_mesh(64)
        return
    end

    if get(ARGS, 1, "") == "sweep"
        arg_p = get(ARGS, 2, "all")
        ps = arg_p == "all" ? SWEEP_PS : (parse(Int, arg_p),)
        outbase = length(ARGS) >= 3 ? ARGS[3] : default_outbase()
        run_sweep(ps, outbase)
        return
    end

    p    = parse(Int, get(ARGS, 1, "3"))
    ndiv = parse(Int, get(ARGS, 2, "10"))

    ensure_mesh(ndiv)
    r = run_case(p, ndiv)

    @printf("Simply supported plate under sinusoidal load (Gauss + Nitsche)\n")
    @printf("  p = %d, ndiv = %d, nodes = %d\n", r.p, r.ndiv, r.nʷ)
    @printf("  h = %.6e, support = %.2f h\n", r.h, support_ratio(r.p))
    @printf("  L2 = %.6e\n", r.L2)
    @printf("  H1 = %.6e\n", r.H1)
    @printf("  H2 = %.6e   (He, relative energy norm of w - wʰ)\n", r.He)
end

main()
