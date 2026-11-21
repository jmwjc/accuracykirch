# ---------------------------------------------------------------------------
# sinusoidal_rkgsi.jl
#
# Simply supported rectangular plate under a sinusoidal body force,
# solved with the mixed Hellinger--Reissner RKGSI formulation.
#
# Exact solution (Navier):
#
#     q̄ = q̄₀ sin(αx) sin(βy),   α = π/a,   β = π/b
#     w = q̄₀ / (D (α² + β²)²) sin(αx) sin(βy)
#
# Boundary conditions:
#
#     w = 0,   Mₙ = 0   on Γ
#
# The deflection w is the RK scalar field; the bending moment Mᵢⱼ is an
# independent piecewise polynomial of degree q - 2 on each integration cell.
# The mixed saddle-point system is assembled and solved directly.
#
# Usage:
#
#     julia program/sinusoidal_rkgsi.jl [p] [q] [ndiv]
#
#     p     RK basis degree: 2, 3 or 4          (default 3)
#     q     piecewise moment degree: 0, 1, 2    (default 2)
#     ndiv  number of subdivisions per edge     (default 10)
#
# The integration consistency order is q + 2, so q = 0, 1, 2 correspond to
# constant, linear and quadratic moment fields, respectively.
#
# The mesh is generated once into program/msh/ by ensure_mesh and is only read
# inside run_case; no mesh is generated during the solve.
#
# The piecewise constant field (q = 0) leaves the domain coupling
# ∫ Mᵢⱼ,ᵢⱼ w dΩ empty.  A penalty term αw ∫ w δw dΓ is added to the w block
# in that case to enforce the simply supported condition and remove the
# singularity of the saddle-point system.
#
# The reported H2 error is the package's relative energy norm of the
# deflection error, He = ‖w - wʰ‖_E / ‖wʰ‖_E.  The recovered-moment energy
# norm is still available internally but is not used in the CSV statistics.
#
# Run from the project root:
#
#     julia program/sinusoidal_rkgsi.jl 3 2 10
# ---------------------------------------------------------------------------

# macOS: Homebrew Julia, Gmsh and UMFPACK may each bundle libomp.
# Allow the duplicate runtime; see OMP Error #15.
ENV["KMP_DUPLICATE_LIB_OK"] = "TRUE"

using ApproxOperator
import ApproxOperator: get𝑛𝑝
import ApproxOperator.GmshImport:
    getPhysicalGroups, get𝑿ᵢ, getElements,
    getPiecewiseElements, getPiecewiseBoundaryElements
import ApproxOperator.KirchhoffPlate:
    ∫MC⁻¹MdΩ, ∫∇²MwdΩ,
    ∫VwdΓ, ∫MθdΓ, wP,
    ∫wqdΩ, L₂, H₁, Hₑ
import ApproxOperator.MindlinPlate: ∫αwwdΓ
import BenchmarkExample
import Gmsh: gmsh

using LinearAlgebra
using SparseArrays
using Printf
using Dates

# Dictionary-backed sparse accumulator.  The weak-form operators assemble with
# `k[i,j] += v`; writing that directly into a SparseMatrixCSC repeatedly grows
# its internal arrays and becomes prohibitively slow on fine meshes.  The COO
# accumulator stores the entries in a Dict and is converted to CSC once.
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
# Physical and geometric parameters (edit as needed)
# ---------------------------------------------------------------------------
const D      = 1.0            # bending stiffness
const ν      = 0.3            # Poisson ratio
const L      = 1.0            # edge length of the square plate
const q̄₀     = 4.0*π^4       # load amplitude, so that w = sin(πx) sin(πy)
const post_order = 12         # quadrature order for error evaluation

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

M₁₁(x, y)    = D*A*(α^2 + ν*β^2)*sin(α*x)*sin(β*y)
M₂₂(x, y)    = D*A*(ν*α^2 + β^2)*sin(α*x)*sin(β*y)
M₁₂(x, y)    = -D*(1.0 - ν)*A*α*β*cos(α*x)*cos(β*y)
q̄(x, y)      = q̄₀*sin(α*x)*sin(β*y)

# 5-point Gauss--Lobatto rule, used for cell boundaries and corners
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
# Helpers
# ---------------------------------------------------------------------------
const mshdir = joinpath(@__DIR__, "msh")

# Path of the mesh used by one run.  The mesh itself is never generated here:
# generate it once with ensure_mesh (or generate all meshes before a sweep).
function meshfile(ndiv::Int)
    tagL = replace(string(L), "." => "p")
    return joinpath(mshdir, "sinusoidal_L$(tagL)_$(ndiv).msh")
end

# Generate a mesh in program/msh/ if it does not exist yet.
function ensure_mesh(ndiv::Int)
    mkpath(mshdir)
    file = meshfile(ndiv)
    if !isfile(file)
        generate_mesh(file, ndiv)
    end
    return file
end

# Square mesh on (0,L) x (0,L) with physical groups Γ¹..Γ⁴, Γ and Ω.
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

function wtype(p::Int)
    return p == 2 ? ReproducingKernel{:Quadratic2D,:□,:QuinticSpline} :
           p == 3 ? ReproducingKernel{:Cubic2D,:□,:QuinticSpline} :
                    ReproducingKernel{:Quartic2D,:□,:QuinticSpline}
end

function mtype(q::Int)
    return q == 0 ? PiecewisePolynomial{:Constant} :
           q == 1 ? PiecewisePolynomial{:Linear2D} :
           q == 2 ? PiecewisePolynomial{:Quadratic2D} :
                    error("q must be 0, 1 or 2 for the piecewise moment field")
end

function support_ratio(p::Int)
    return p == 2 ? 2.5 : p == 3 ? 3.5 : 4.5
end

# Assemble the H2-equivalent energy norm of the recovered moment error.
function moment_error(elements_m, d, nM)
    C₁₁ = 1.0/(D*(1.0 - ν^2))
    C₁₂ = -ν/(D*(1.0 - ν^2))
    C₃₃ = 2.0/(D*(1.0 - ν))

    ΔW² = 0.0
    W̄²  = 0.0
    for (e, elem) in enumerate(elements_m)
        off = 3*(e - 1)*nM
        for ξ in elem.𝓖
            N = ξ[:𝝭]
            m₁₁ = 0.0; m₂₂ = 0.0; m₁₂ = 0.0
            for i in 1:nM
                m₁₁ += N[i]*d[off + 3i - 2]
                m₂₂ += N[i]*d[off + 3i - 1]
                m₁₂ += N[i]*d[off + 3i]
            end
            M̄₁₁ = M₁₁(ξ.x, ξ.y)
            M̄₂₂ = M₂₂(ξ.x, ξ.y)
            M̄₁₂ = M₁₂(ξ.x, ξ.y)
            Δ₁₁ = m₁₁ - M̄₁₁
            Δ₂₂ = m₂₂ - M̄₂₂
            Δ₁₂ = m₁₂ - M̄₁₂
            ΔW² += (C₁₁*Δ₁₁^2 + 2.0*C₁₂*Δ₁₁*Δ₂₂ +
                    C₁₁*Δ₂₂^2 + C₃₃*Δ₁₂^2)*ξ.𝑤
            W̄²  += (C₁₁*M̄₁₁^2 + 2.0*C₁₂*M̄₁₁*M̄₂₂ +
                    C₁₁*M̄₂₂^2 + C₃₃*M̄₁₂^2)*ξ.𝑤
        end
    end
    return sqrt(ΔW²/W̄²), ΔW², W̄²
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
    h = 1.0/ndiv
    support = support_ratio(p)*h
    αw = 1.0e3 * D / h^3
    push!(nodes, :s₁ => fill(support, nʷ),
                 :s₂ => fill(support, nʷ),
                 :s₃ => fill(support, nʷ))

    sp = RegularGrid(x, y, z; n = 3, γ = 5)
    order = p == 2 ? 4 : p == 3 ? 4 : 6

    # --- w field -----------------------------------------------------------
    elements_w = getElements(nodes, entities["Ω"], wtype(p), order, sp)
    prescribe!(elements_w, :q => q̄)
    set∇²𝝭!(elements_w)

    # --- independent moment field -----------------------------------------
    elements_m = getPiecewiseElements(entities["Ω"], mtype(q), order)
    if q <= 1
        set∇𝝭!(elements_m)
        push!(elements_m, :∂²𝝭∂x², :∂²𝝭∂y², :∂²𝝭∂x∂y)
    else
        set∇²𝝭!(elements_m)
    end

    nₑ = length(elements_m)
    nM = get𝑛𝑝(elements_m[1])
    nᵐ = nₑ*nM

    kᵐᵐ = COO(3*nᵐ, 3*nᵐ)
    kᵐʷ = COO(3*nᵐ, nʷ)
    fʷ  = zeros(nʷ)
    fᵐ  = zeros(3*nᵐ)

    (∫MC⁻¹MdΩ => elements_m)(kᵐᵐ; D = D, ν = ν)
    (∫∇²MwdΩ  => (elements_m, elements_w))(kᵐʷ)
    (∫wqdΩ    => elements_w)(fʷ)

    # --- cell-boundary terms Γ_C ------------------------------------------
    elements_Γ_w = getElements(nodes, entities["Γ"], wtype(p), GL5, sp;
                               normal = true)
    elements_Γ_m = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"],
                                                mtype(q), GL5)
    set∇𝝭!(elements_Γ_w)
    set∇𝝭!(elements_Γ_m)
    for op in (∫VwdΓ  => (elements_Γ_m, elements_Γ_w),
               ∫MθdΓ => (elements_Γ_m, elements_Γ_w),
               wP    => (elements_Γ_m, elements_Γ_w))
        op(kᵐʷ)
    end

    # --- essential condition w = 0 on Γ ----------------------------------
    elements_Γ¹_w = getElements(nodes, entities["Γ¹"], wtype(p), GL5, sp; normal = true)
    elements_Γ²_w = getElements(nodes, entities["Γ²"], wtype(p), GL5, sp; normal = true)
    elements_Γ³_w = getElements(nodes, entities["Γ³"], wtype(p), GL5, sp; normal = true)
    elements_Γ⁴_w = getElements(nodes, entities["Γ⁴"], wtype(p), GL5, sp; normal = true)
    elements_Γ¹_m = getElements(entities["Γ¹"], entities["Γ"], elements_Γ_m)
    elements_Γ²_m = getElements(entities["Γ²"], entities["Γ"], elements_Γ_m)
    elements_Γ³_m = getElements(entities["Γ³"], entities["Γ"], elements_Γ_m)
    elements_Γ⁴_m = getElements(entities["Γ⁴"], entities["Γ"], elements_Γ_m)
    for e in (elements_Γ¹_w, elements_Γ²_w, elements_Γ³_w, elements_Γ⁴_w)
        prescribe!(e, :w => w, :g => w)
        set∇𝝭!(e)
    end
    elements_Γᵍ_w = elements_Γ¹_w ∪ elements_Γ²_w ∪ elements_Γ³_w ∪ elements_Γ⁴_w
    elements_Γᵍ_m = elements_Γ¹_m ∪ elements_Γ²_m ∪ elements_Γ³_m ∪ elements_Γ⁴_m
    (∫VwdΓ => (elements_Γᵍ_m, elements_Γᵍ_w))(kᵐʷ, fᵐ)
    (wP    => (elements_Γᵍ_m, elements_Γᵍ_w))(kᵐʷ, fᵐ)

    # --- solve the saddle-point system ------------------------------------
    # Penalty stabilization of the w block for the constant moment field.
    kʷʷ = COO(nʷ, nʷ)
    if q == 0
        (∫αwwdΓ => elements_Γᵍ_w)(kʷʷ, fʷ; α = αw)
    end

    kmm = to_sparse(kᵐᵐ)
    kmw = to_sparse(kᵐʷ)
    kww = to_sparse(kʷʷ)
    K = [kmm kmw; kmw' kww]
    d = K \ [fᵐ; fʷ]
    w_sol = d[3*nᵐ+1:3*nᵐ+nʷ]

    # --- H2-equivalent error from the recovered moment --------------------
    elements_m_post = getPiecewiseElements(entities["Ω"], mtype(q), post_order)
    set𝝭!(elements_m_post)
    H2, _, _ = moment_error(elements_m_post, d, nM)

    # --- L2, H1 and deflection energy errors ------------------------------
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
            L2, H1, H2, He)
end

# ---------------------------------------------------------------------------
# Sweep and CSV output
# ---------------------------------------------------------------------------
const SWEEP_PS    = (2, 3, 4)
const SWEEP_QS    = (0, 1, 2)
const SWEEP_NDIVS = (8, 16, 32, 64)

function timestamp()
    return Dates.format(now(), "yyyymmdd_HHMMSS")
end

function default_outbase()
    return joinpath(@__DIR__, "results", "sinusoidal_rkgsi_" * timestamp())
end

function ensure_new_output(outbase)
    for suffix in ("", "_rates", "_config")
        file = outbase * suffix * ".csv"
        if isfile(file)
            error("refusing to overwrite existing result: $file; choose a new output prefix")
        end
    end
end

function run_sweep(ps = SWEEP_PS, qs = SWEEP_QS,
                   outbase = default_outbase())
    # Generate every mesh once, before the solve loop.
    ensure_mesh(8)
    ensure_mesh(16)
    ensure_mesh(32)
    ensure_mesh(64)

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
            push!(rows, (; p, q, ndiv, nodes = 0, h = "", L2 = "",
                            H1 = "", H2 = "", He = "",
                            valid = "no", note = "run failed"))
            continue
        end

        valid = isfinite(r.L2) && isfinite(r.H1) && isfinite(r.He) &&
                r.L2 < 1.0 && r.H1 < 1.0 && r.He < 1.0
        note = valid ? "" : "singular / invalid"

        push!(rows, (; p, q, ndiv,
                        nodes = r.nʷ, h = r.h,
                        L2 = r.L2, H1 = r.H1, H2 = r.He, He = r.He,
                        valid = valid ? "yes" : "no", note))
    end

    write_results(rows, outbase)
    return rows
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
        n = length(cols[1])
        for i in 1:n
            println(io, join([csv_cell(cols[j][i]) for j in eachindex(cols)], ","))
        end
    end
    return file
end

function write_results(rows, outbase = default_outbase())
    ensure_new_output(outbase)

    raw = (
        p      = Int[r.p for r in rows],
        q      = Int[r.q for r in rows],
        q_int  = Int[r.q + 2 for r in rows],
        ndiv   = Int[r.ndiv for r in rows],
        nodes  = Int[r.nodes for r in rows],
        h      = Any[r.h for r in rows],
        L2     = Any[r.L2 for r in rows],
        H1     = Any[r.H1 for r in rows],
        H2     = Any[r.H2 for r in rows],
        He     = Any[r.He for r in rows],
        valid  = String[r.valid for r in rows],
        note   = String[r.note for r in rows],
    )

    rate_p = Int[]
    rate_q = Int[]
    norm   = String[]
    r_8_16 = Any[]
    r_16_32 = Any[]
    r_32_64 = Any[]
    r_all  = Any[]

    for p in unique(r.p for r in rows), q in unique(r.q for r in rows)
        sub = sort([r for r in rows if r.p == p && r.q == q], by = r -> r.ndiv)
        length(sub) == 4 || continue
        ok = all(r -> r.valid == "yes", sub)

        for (name, key) in (("L2", :L2), ("H1", :H1), ("H2", :H2))
            vals = [getproperty(r, key) for r in sub]
            push!(rate_p, p)
            push!(rate_q, q)
            push!(norm, name)
            if ok && all(isfinite, vals) && all(v -> v > 0, vals)
                push!(r_8_16,  log(vals[1]/vals[2])/log(2))
                push!(r_16_32, log(vals[2]/vals[3])/log(2))
                push!(r_32_64, log(vals[3]/vals[4])/log(2))
                push!(r_all,   log(vals[1]/vals[4])/log(8))
            else
                push!(r_8_16, "")
                push!(r_16_32, "")
                push!(r_32_64, "")
                push!(r_all, "")
            end
        end
    end

    rates = (;
        p = rate_p,
        q = rate_q,
        norm,
        rate_8_16 = r_8_16,
        rate_16_32 = r_16_32,
        rate_32_64 = r_32_64,
        rate_all = r_all,
    )

    config = (
        item = ["D", "nu", "L", "q0", "q definition", "error norm", "support"],
        value = [string(D), string(ν), string(L), string(q̄₀),
                 "q = piecewise moment degree; integration order = q + 2",
                 "H2 = He = relative energy norm of w - w^h",
                 "p=2: 2.5h; p=3: 3.5h; p=4: 4.5h"],
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
        arg_q = get(ARGS, 3, "all")
        ps = arg_p == "all" ? SWEEP_PS : (parse(Int, arg_p),)
        qs = arg_q == "all" ? SWEEP_QS : (parse(Int, arg_q),)
        outbase = length(ARGS) >= 4 ? ARGS[4] : default_outbase()
        run_sweep(ps, qs, outbase)
        return
    end

    p    = parse(Int, get(ARGS, 1, "3"))
    q    = parse(Int, get(ARGS, 2, "2"))
    ndiv = parse(Int, get(ARGS, 3, "10"))

    # Mesh generation is kept outside the solve, one file per ndiv in program/msh/.
    ensure_mesh(ndiv)

    r = run_case(p, q, ndiv)

    @printf("Simply supported plate under sinusoidal load (HR-RKGSI)\n")
    @printf("  p = %d, q = %d, ndiv = %d\n", r.p, r.q, r.ndiv)
    @printf("  nodes = %d, cells = %d, M basis = %d, M dofs = %d\n",
            r.nʷ, r.nₑ, r.nM, r.n_moment_dofs)
    @printf("  h = %.6e, support = %.2f h\n", r.h, support_ratio(r.p))
    @printf("  L2 = %.6e\n", r.L2)
    @printf("  H1 = %.6e\n", r.H1)
    @printf("  H2 = %.6e   (He, relative energy norm of w - wʰ)\n", r.He)
end

main()
