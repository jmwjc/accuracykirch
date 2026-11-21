# Merge the per-p XLSX files produced by
#     julia program/sinusoidal.jl sweep <p> all program/results/p<p>_support.xlsx
# into a single workbook.
using XLSX

files = isempty(ARGS) ? ["program/results/p2_support.xlsx",
                         "program/results/p3_support.xlsx",
                         "program/results/p4_support.xlsx"] : ARGS
outfile = "program/sinusoidal_results.xlsx"

function read_cols(file, sheet)
    t = XLSX.readtable(file, sheet)
    return t.column_labels, [Any[x for x in col] for col in t.data]
end

function append_cols!(dst, src)
    for i in eachindex(dst)
        append!(dst[i], src[i])
    end
end

function nt(labels, cols)
    return NamedTuple{Tuple(labels)}(Tuple(cols))
end

function dedup_config!(labels, cols)
    i = findfirst(==(:item), labels)
    i === nothing && return
    seen = Set{Any}()
    keep = Bool[]
    for item in cols[i]
        if item in seen
            push!(keep, false)
        else
            push!(seen, item)
            push!(keep, true)
        end
    end
    for j in eachindex(cols)
        cols[j] = cols[j][keep]
    end
end

function sort_rows!(labels, cols, keys, normorder = nothing)
    idx = [findfirst(==(k), labels) for k in keys]
    n = length(cols[1])
    order = sortperm(1:n, by = i -> (
        cols[idx[1]][i], cols[idx[2]][i],
        normorder === nothing ? cols[idx[3]][i] :
            findfirst(==(cols[idx[3]][i]), normorder)))
    for col in cols
        col .= col[order]
    end
end

labels = nothing
cols = nothing
for f in files
    isfile(f) || continue
    l, c = read_cols(f, "raw")
    labels === nothing && (labels = l; cols = [Any[] for _ in l])
    append_cols!(cols, c)
end

# Blank error columns of invalid rows.
iv = findfirst(==(:valid), labels)
if iv !== nothing
    for j in eachindex(cols[1])
        if cols[iv][j] == "no"
            for name in (:L2, :H1, :H2, :He)
                k = findfirst(==(name), labels)
                k === nothing || (cols[k][j] = "")
            end
        end
    end
end
sort_rows!(labels, cols, (:p, :q, :ndiv))
raw = nt(labels, cols)

labels = nothing
cols = nothing
for f in files
    isfile(f) || continue
    l, c = read_cols(f, "rates")
    labels === nothing && (labels = l; cols = [Any[] for _ in l])
    append_cols!(cols, c)
end
sort_rows!(labels, cols, (:p, :q, :norm), ["L2", "H1", "H2"])
rates = nt(labels, cols)

labels = nothing
cols = nothing
for f in files
    isfile(f) || continue
    l, c = read_cols(f, "config")
    labels === nothing && (labels = l; cols = [Any[] for _ in l])
    append_cols!(cols, c)
end
dedup_config!(labels, cols)
config = nt(labels, cols)

XLSX.writetable(outfile,
                "raw" => raw,
                "rates" => rates,
                "config" => config;
                overwrite = true)
println("written: ", outfile)
