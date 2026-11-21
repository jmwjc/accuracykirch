# ---------------------------------------------------------------------------
# circle_rkgsi.jl
#
# Mixed Hellinger--Reissner RKGSI test for the circular plate described in
# manuscript/sections/circle.tex.
#
# The computation uses the full circular plate.  The full disk is generated
# as a single surface, so no artificial symmetry edges or center corner are
# present in the mesh.
#
# Exact solution in circle.tex:
#
#     w = q̄ (R² - r²) / (64 D) * ((5 + ν)/(1 + ν) R² - r²)
#
# This closed form has a non-zero slope at r = R, so it is the
# simply-supported solution used in circle.tex.
# Outer circle Γᵒ: w = 0, Mₙ = M̄ (natural).
#
# Run from the project root:
#
#     julia program/circle_rkgsi.jl [p] [q] [ndiv]
#     julia program/circle_rkgsi.jl sweep [p] [q] [outbase]
# ---------------------------------------------------------------------------

ENV["KMP_DUPLICATE_LIB_OK"] = "TRUE"

using ApproxOperator
import ApproxOperator: get𝑛𝑝
import ApproxOperator.GmshImport:
    getPhysicalGroups, get𝑿ᵢ, getElements,
    getPiecewiseElements, getPiecewiseBoundaryElements
import ApproxOperator.KirchhoffPlate:
    ∫MC⁻¹MdΩ, ∫∇²MwdΩ,
    ∫VwdΓ, ∫MθdΓ, wP,
    ∫θₙMₙₙdΓ,
    ∫wqdΩ, L₂, H₁, Hₑ
import ApproxOperator.MindlinPlate: ∫αwwdΓ
import BenchmarkExample
import Gmsh: gmsh

using LinearAlgebra
using SparseArrays
using Printf
using Dates

# Dictionary-backed sparse accumulator, as in sinusoidal_rkgsi.jl.
struct COO <: AbstractMatrix{Float64}
    nrow::Int
    ncol::Int
    vals::Dict{Tuple{Int,Int},Float64}
end
COO(nrow::Int, ncol::Int) = COO(nrow, ncol, Dict{Tuple{Int,Int},Float64}())
Base.size(A::COO) = (A.nrow, A.ncol)
Base.getindex(A::COO, i::Int, j::Int) = get(A.vals, (i, j), 0.0)
Base.setindex!(A::COO, v::Float64, i::Int, j::Int) = (A.vals[(i, j)] = v)

function to_sparse(A::COO)
    n = length(A.vals)
    I = Vector{Int}(undef, n)
    J = Vector{Int}(undef, n)
    V = Vector{Float64}(undef, n)
    k = 0
    for ((i, j), v) in A.vals
        k += 1
        I[k] = i
        J[k] = j
        V[k] = v
    end
    return sparse(I, J, V, A.nrow, A.ncol)
end

# ---------------------------------------------------------------------------
# Problem parameters
# ---------------------------------------------------------------------------
const R       = 5.0
const D       = 2.0e4
const ν       = 0.3
const q̄₀     = 1.28e6
const post_order = 12

const Acoef = (5.0 + ν)/(1.0 + ν)
const Ccoef = q̄₀/(64.0*D)

r²(x, y) = x^2 + y^2

# Simply supported closed form (circle.tex)
w(x, y)      = Ccoef*(R^2 - r²(x, y))*(Acoef*R^2 - r²(x, y))
w₁(x, y)     = 2.0*Ccoef*x*(-(1.0 + Acoef)*R^2 + 2.0*r²(x, y))
w₂(x, y)     = 2.0*Ccoef*y*(-(1.0 + Acoef)*R^2 + 2.0*r²(x, y))
w₁₁(x, y)    = 2.0*Ccoef*(-(1.0 + Acoef)*R^2 + 2.0*r²(x, y) + 4.0*x^2)
w₁₂(x, y)    = 8.0*Ccoef*x*y
w₂₂(x, y)    = 2.0*Ccoef*(-(1.0 + Acoef)*R^2 + 2.0*r²(x, y) + 4.0*y^2)

# Exact bending moments of the selected closed form.
M₁₁(x, y) = -D*(w₁₁(x, y) + ν*w₂₂(x, y))
M₂₂(x, y) = -D*(ν*w₁₁(x, y) + w₂₂(x, y))
M₁₂(x, y) = -D*(1.0 - ν)*w₁₂(x, y)

# Normal bending moment, used for the natural boundary value M̄.
M(x, y, z, n₁, n₂) = M₁₁(x, y)*n₁^2 + M₂₂(x, y)*n₂^2 +
                     2.0*M₁₂(x, y)*n₁*n₂

q̄(x, y) = q̄₀
θ(x, y, n₁, n₂) = w₁(x, y)*n₁ + w₂(x, y)*n₂

# 5-point Gauss-Lobatto, including endpoints for corner terms.
const wq = [1.0/10.0, 49.0/90.0, 32.0/45.0, 49.0/90.0, 1.0/10.0]
const GL5 = (Float64[
    -1.0, 0.0, 0.0,
    -sqrt(3.0/7.0), 0.0, 0.0,
     0.0, 0.0, 0.0,
     sqrt(3.0/7.0), 0.0, 0.0,
     1.0, 0.0, 0.0,
], wq)

# ---------------------------------------------------------------------------
# Mesh: full disk
# ---------------------------------------------------------------------------
const mshdir = joinpath(@__DIR__, "msh")

function meshfile(ndiv::Int)
    return joinpath(mshdir, "circle_full_ref_$(ndiv).msh")
end

const BASE_NDIV = 4

function refine_level(ndiv::Int)
    n = round(Int, log2(ndiv/BASE_NDIV))
    ndiv == BASE_NDIV*2^n || error("ndiv must be BASE_NDIV * 2^n")
    0 <= n <= 3 || error("only 0 to 3 refinements are supported")
    return n
end

# Quarter-disk mesh with a single surface.  The geometry follows
# BenchmarkExample.Circular, but the original three-surface patch is merged
# into one surface so that ApproxOperator sees a single physical surface.
function generate_circle_mesh(file::String, ndiv::Int)
    nref = refine_level(ndiv)
    gmsh.initialize()
    gmsh.model.add("Circle full")

    # Full disk: always start from the same coarse mesh and refine it.
    lc = R/BASE_NDIV
    p₀ = gmsh.model.geo.addPoint(0.0, 0.0, 0.0, lc, 1)
    p₁ = gmsh.model.geo.addPoint(R, 0.0, 0.0, lc, 2)
    p₂ = gmsh.model.geo.addPoint(0.0, R, 0.0, lc, 3)
    p₃ = gmsh.model.geo.addPoint(-R, 0.0, 0.0, lc, 4)
    p₄ = gmsh.model.geo.addPoint(0.0, -R, 0.0, lc, 5)

    Γ₁ = gmsh.model.geo.addCircleArc(2, 1, 3, 1)
    Γ₂ = gmsh.model.geo.addCircleArc(3, 1, 4, 2)
    Γ₃ = gmsh.model.geo.addCircleArc(4, 1, 5, 3)
    Γ₄ = gmsh.model.geo.addCircleArc(5, 1, 2, 4)

    gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)
    Ω = gmsh.model.geo.addPlaneSurface([1], 1)
    gmsh.model.geo.synchronize()

    gmsh.model.addPhysicalGroup(1, [Γ₁, Γ₂, Γ₃, Γ₄], -1, "Γᵒ")
    gmsh.model.addPhysicalGroup(2, [Ω], -1, "Ω")

    gmsh.model.mesh.setAlgorithm(2, 1, 1)
    gmsh.model.mesh.generate(2)
    for _ in 1:nref
        gmsh.model.mesh.refine()
    end
    gmsh.model.mesh.setOrder(1)

    tag = BenchmarkExample.addEdgeElements((2, 1), 1)
    gmsh.model.geo.addPhysicalGroup(1, [tag], -1, "Γ")
    gmsh.model.geo.synchronize()
    gmsh.write(file)
    gmsh.finalize()
end

function ensure_mesh(ndiv::Int)
    mkpath(mshdir)
    file = meshfile(ndiv)
    if !isfile(file)
        generate_circle_mesh(file, ndiv)
    end
    return file
end

function wtype(p::Int)
    return p == 2 ? ReproducingKernel{:Quadratic2D,:□,:QuinticSpline} :
           p == 3 ? ReproducingKernel{:Cubic2D,:□,:QuinticSpline} :
                    ReproducingKernel{:Quartic2D,:□,:QuinticSpline}
end

function mtype(q::Int)
    return q == 0 ? PiecewisePolynomial{:Constant} :
           q == 1 ? PiecewisePolynomial{:Linear2D} :
           q == 2 ? PiecewisePolynomial{:Quadratic2D} :
                    error("q must be 0, 1 or 2")
end

function support_ratio(p::Int)
    return p == 2 ? 2.5 : p == 3 ? 3.5 : 4.5
end

# ---------------------------------------------------------------------------
# One run
# ---------------------------------------------------------------------------
function run_case(p::Int, q::Int, ndiv::Int)
    p in (2, 3, 4) || error("p must be 2, 3 or 4")
    q in (0, 1, 2) || error("q must be 0, 1 or 2")
    ndiv > 1       || error("ndiv must be larger than one")

    file = meshfile(ndiv)
    isfile(file) || error("mesh not found: $file; call ensure_mesh($ndiv) first")

    gmsh.initialize()
    gmsh.open(file)
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()

    x, y, z = nodes.x, nodes.y, nodes.z
    nʷ = length(nodes)

    sp = RegularGrid(x, y, z; n = 3, γ = 5)
    order = p == 2 ? 4 : p == 3 ? 4 : 6

    # Number of triangular cells is needed for the equivalent element size
    # of the acute-triangle mesh.
    elements_m = getPiecewiseElements(entities["Ω"], mtype(q), order)
    nₑ = length(elements_m)
    A_total = π*R^2
    h = sqrt(4.0*A_total/(sqrt(3.0)*nₑ))

    support = support_ratio(p)*h
    αw = 1.0e3 * D / h^3
    push!(nodes, :s₁ => fill(support, nʷ),
                 :s₂ => fill(support, nʷ),
                 :s₃ => fill(support, nʷ))

    # --- domain -----------------------------------------------------------
    elements_w = getElements(nodes, entities["Ω"], wtype(p), order, sp)
    prescribe!(elements_w, :q => q̄)
    set∇²𝝭!(elements_w)

    if q <= 1
        set∇𝝭!(elements_m)
        push!(elements_m, :∂²𝝭∂x², :∂²𝝭∂y², :∂²𝝭∂x∂y)
    else
        set∇²𝝭!(elements_m)
    end

    nM = get𝑛𝑝(elements_m[1])
    nᵐ = nₑ*nM

    kᵐᵐ = COO(3*nᵐ, 3*nᵐ)
    kᵐʷ = COO(3*nᵐ, nʷ)
    fʷ  = zeros(nʷ)
    fᵐ  = zeros(3*nᵐ)

    (∫MC⁻¹MdΩ => elements_m)(kᵐᵐ; D = D, ν = ν)
    (∫∇²MwdΩ  => (elements_m, elements_w))(kᵐʷ)
    (∫wqdΩ    => elements_w)(fʷ)

    # --- cell-boundary Γ_C terms -----------------------------------------
    elements_Γ_w = getElements(nodes, entities["Γ"], wtype(p), GL5, sp;
                               normal = true)
    elements_Γ_m = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"],
                                                mtype(q), GL5)
    set∇𝝭!(elements_Γ_w)
    q <= 1 ? set∇𝝭!(elements_Γ_m) : set∇²𝝭!(elements_Γ_m)
    for op in (∫VwdΓ  => (elements_Γ_m, elements_Γ_w),
               ∫MθdΓ => (elements_Γ_m, elements_Γ_w),
               wP    => (elements_Γ_m, elements_Γ_w))
        op(kᵐʷ)
    end

    # --- essential boundary: simply supported outer circle ---------------
    # Γ_w: w = 0 and natural Mₙ = M̄.
    Γ_w_elems = Any[]

    for tag in entities["Γᵒ"].second
        Γw = getElements(nodes, 1 => [tag], wtype(p), GL5, sp; normal = true)
        Γm = getElements(1 => [tag], entities["Γ"], elements_Γ_m)
        prescribe!(Γw, :w => w, :g => w, :M => M)
        set∇𝝭!(Γw)
        (∫VwdΓ => (Γm, Γw))(kᵐʷ, fᵐ)
        (wP    => (Γm, Γw))(kᵐʷ, fᵐ)
        push!(Γ_w_elems, Γw)

        # Natural boundary Mₙ = M̄ on the outer circle.
        (∫θₙMₙₙdΓ => Γw)(fʷ)
    end

    # --- penalty for the constant moment field (q = 0) -------------------
    kʷʷ = COO(nʷ, nʷ)
    if q == 0
        for Γw in Γ_w_elems
            (∫αwwdΓ => Γw)(kʷʷ, fʷ; α = αw)
        end
    end

    kmm = to_sparse(kᵐᵐ)
    kmw = to_sparse(kᵐʷ)
    kww = to_sparse(kʷʷ)
    K = [kmm kmw; kmw' kww]
    d = K \ [fᵐ; fʷ]
    w_sol = d[3*nᵐ+1:3*nᵐ+nʷ]

    # --- errors -----------------------------------------------------------
    push!(nodes, :d => w_sol)
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

    return (; p, q, ndiv, h, nʷ, nₑ, nM,
            n_moment_dofs = 3*nᵐ,
            L2, H1, H2 = He, He)
end

# ---------------------------------------------------------------------------
# Sweep / CSV output (timestamped, refuses to overwrite)
# ---------------------------------------------------------------------------
const SWEEP_PS    = (2, 3, 4)
const SWEEP_QS    = (0, 1, 2)
# Four nested meshes: the same base mesh, then 0, 1, 2 and 3 Gmsh
# refinements.  The effective ndiv values are 4, 8, 16 and 32.
const SWEEP_NDIVS = (BASE_NDIV*2^0, BASE_NDIV*2^1,
                     BASE_NDIV*2^2, BASE_NDIV*2^3)

function timestamp()
    return Dates.format(now(), "yyyymmdd_HHMMSS")
end

function default_outbase()
    return joinpath(@__DIR__, "results", "circle_rkgsi_" * timestamp())
end

function ensure_new_output(outbase)
    for suffix in ("", "_rates", "_config")
        f = outbase * suffix * ".csv"
        isfile(f) && error("refusing to overwrite existing result: $f")
    end
end

function csv_cell(x)
    s = ismissing(x) ? "" : string(x)
    (occursin(',', s) || occursin('"', s) || occursin('\n', s)) &&
        return "\"" * replace(s, "\"" => "\"\"") * "\""
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
        p = Int[r.p for r in rows], q = Int[r.q for r in rows],
        q_int = Int[r.q + 2 for r in rows],
        ndiv = Int[r.ndiv for r in rows],
        nodes = Int[r.nʷ for r in rows],
        h = Any[r.h for r in rows],
        L2 = Any[r.L2 for r in rows], H1 = Any[r.H1 for r in rows],
        H2 = Any[r.H2 for r in rows], He = Any[r.He for r in rows],
        valid = String[isfinite(r.L2) && isfinite(r.H1) && isfinite(r.He) ? "yes" : "no"
                       for r in rows],
    )

    rate_p = Int[]; rate_q = Int[]; rate_norm = String[]
    r1 = Any[]; r2 = Any[]; r3 = Any[]; ra = Any[]
    for p in unique(r.p for r in rows), q in unique(r.q for r in rows)
        sub = sort([r for r in rows if r.p == p && r.q == q], by = r -> r.ndiv)
        length(sub) == 4 || continue
        ok = all(r -> r.L2 isa Real && r.H1 isa Real && r.He isa Real &&
                       r.h isa Real && r.h > 0 &&
                       isfinite(r.L2) && isfinite(r.H1) && isfinite(r.He) &&
                       r.L2 > 0 && r.H1 > 0 && r.He > 0, sub)
        for (name, key) in (("L2", :L2), ("H1", :H1), ("H2", :H2))
            push!(rate_p, p); push!(rate_q, q); push!(rate_norm, name)
            if ok
                vals = [getproperty(r, key) for r in sub]
                hs = [r.h for r in sub]
                push!(r1, log(vals[1]/vals[2])/log(hs[1]/hs[2]))
                push!(r2, log(vals[2]/vals[3])/log(hs[2]/hs[3]))
                push!(r3, log(vals[3]/vals[4])/log(hs[3]/hs[4]))
                push!(ra, log(vals[1]/vals[4])/log(hs[1]/hs[4]))
            else
                for v in (r1, r2, r3, ra); push!(v, ""); end
            end
        end
    end

    rates = (; p = rate_p, q = rate_q, norm = rate_norm,
               rate_4_8 = r1, rate_8_16 = r2, rate_16_32 = r3, rate_all = ra)

    config = (
        item = ["method", "BC", "R", "D", "nu", "q0", "q definition",
                "error norm", "support"],
        value = ["RKGSI (mixed HR)",
                 "Γᵒ: w=0, Mₙ=M̄ (simply supported)",
                 string(R), string(D), string(ν), string(q̄₀),
                 "q = piecewise moment degree; integration order = q + 2",
                 "H2 = He = relative energy norm of w - w^h",
                 "p=2: 2.5h; p=3: 3.5h; p=4: 4.5h"],
    )

    write_csv(outbase * ".csv", String.(keys(raw)),
              [Any[x for x in getproperty(raw, k)] for k in keys(raw)])
    write_csv(outbase * "_rates.csv", String.(keys(rates)),
              [Any[x for x in getproperty(rates, k)] for k in keys(rates)])
    write_csv(outbase * "_config.csv", String.(keys(config)),
              [Any[x for x in getproperty(config, k)] for k in keys(config)])
    @printf("written: %s.csv\n", outbase)
    return outbase
end

function run_sweep(ps = SWEEP_PS, qs = SWEEP_QS, outbase = default_outbase())
    for n in SWEEP_NDIVS; ensure_mesh(n); end
    rows = NamedTuple[]
    for p in ps, q in qs, ndiv in SWEEP_NDIVS
        @printf("running p = %d, q = %d, ndiv = %d ...\n", p, q, ndiv)
        flush(stdout)
        r = try
            run_case(p, q, ndiv)
        catch err
            @warn "run failed" p q ndiv exception = err
            nothing
        end
        if r === nothing
            push!(rows, (; p, q, ndiv, h = "", nʷ = 0,
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
        for n in SWEEP_NDIVS; ensure_mesh(n); end
        return
    end

    if get(ARGS, 1, "") == "sweep"
        arg_p = get(ARGS, 2, "all")
        arg_q = get(ARGS, 3, "all")
        ps = arg_p == "all" ? SWEEP_PS : (parse(Int, arg_p),)
        qs = arg_q == "all" ? SWEEP_QS : (parse(Int, arg_q),)
        outbase = length(ARGS) >= 4 ? ARGS[4] : default_outbase()
        run_sweep(ps, qs, outbase)
        return
    end

    p = parse(Int, get(ARGS, 1, "3"))
    q = parse(Int, get(ARGS, 2, "2"))
    ndiv = parse(Int, get(ARGS, 3, "8"))
    ensure_mesh(ndiv)
    r = run_case(p, q, ndiv)
    @printf("Full circular plate (simply supported): p = %d, q = %d, ndiv = %d\n",
            r.p, r.q, r.ndiv)
    @printf("  nodes = %d, cells = %d, M dofs = %d\n", r.nʷ, r.nₑ, r.n_moment_dofs)
    @printf("  L2 = %.6e\n", r.L2)
    @printf("  H1 = %.6e\n", r.H1)
    @printf("  H2 = %.6e   (He)\n", r.He)
end

main()
