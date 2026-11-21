# Replace the H2 statistic by He in the existing CSV tables and recompute rates.
using Printf

const RAWFILE = "program/sinusoidal_results.csv"
const RATEFILE = "program/sinusoidal_results_rates.csv"
const CONFIGFILE = "program/sinusoidal_results_config.csv"

function parse_cell(s)
    s == "" && return ""
    i = tryparse(Int, s)
    i !== nothing && return i
    v = tryparse(Float64, s)
    return v === nothing ? s : v
end

function read_csv_table(path)
    lines = readlines(path)
    labels = Symbol.(split(lines[1], ','))
    cols = [Any[] for _ in labels]
    for line in lines[2:end]
        isempty(line) && continue
        vals = split(line, ',')
        for j in eachindex(cols)
            push!(cols[j], parse_cell(vals[j]))
        end
    end
    return labels, cols
end

function csv_cell(x)
    s = ismissing(x) ? "" : string(x)
    (occursin(',', s) || occursin('"', s) || occursin('\n', s)) &&
        return "\"" * replace(s, "\"" => "\"\"") * "\""
    return s
end

function write_csv(path, labels, cols)
    open(path, "w") do io
        println(io, join(string.(labels), ","))
        for i in eachindex(cols[1])
            println(io, join([csv_cell(cols[j][i]) for j in eachindex(cols)], ","))
        end
    end
end

labels, cols = read_csv_table(RAWFILE)
iH2 = findfirst(==(:H2), labels)
iHe = findfirst(==(:He), labels)
for i in eachindex(cols[iH2])
    cols[iH2][i] = cols[iHe][i]
end

# Recompute rates for L2, H1 and the new H2 (= He).
ip = findfirst(==(:p), labels)
iq = findfirst(==(:q), labels)
in = findfirst(==(:ndiv), labels)
irate = Dict(:L2 => 1, :H1 => 2, :H2 => 3)
rate_labels = [:p, :q, :norm, :rate_8_16, :rate_16_32, :rate_32_64, :rate_all]
rate_cols = [Any[] for _ in rate_labels]
for p in sort(unique(cols[ip]))
    for q in sort(unique(cols[iq]))
        idx = [i for i in eachindex(cols[ip]) if cols[ip][i] == p && cols[iq][i] == q]
        sort!(idx, by = i -> cols[in][i])
        length(idx) == 4 || continue
        for (name, label) in ((:L2, :L2), (:H1, :H1), (:H2, :H2))
            j = findfirst(==(label), labels)
            vals = [cols[j][i] for i in idx]
            ok = all(v -> v isa Real && isfinite(v) && v > 0, vals)
            push!(rate_cols[1], p)
            push!(rate_cols[2], q)
            push!(rate_cols[3], string(name))
            if ok
                push!(rate_cols[4], log(vals[1]/vals[2])/log(2))
                push!(rate_cols[5], log(vals[2]/vals[3])/log(2))
                push!(rate_cols[6], log(vals[3]/vals[4])/log(2))
                push!(rate_cols[7], log(vals[1]/vals[4])/log(8))
            else
                for k in 4:7; push!(rate_cols[k], ""); end
            end
        end
    end
end

write_csv(RAWFILE, labels, cols)
write_csv(RATEFILE, rate_labels, rate_cols)

# Update the error-norm note in the config table.
clabels, ccols = read_csv_table(CONFIGFILE)
ii = findfirst(==(:item), clabels)
iv = findfirst(==(:value), clabels)
for k in eachindex(ccols[ii])
    if ccols[ii][k] == "error norm"
        ccols[iv][k] = "H2 = He = relative energy norm of w - w^h"
    end
end
write_csv(CONFIGFILE, clabels, ccols)

@printf("updated: %s\n", RAWFILE)
@printf("updated: %s\n", RATEFILE)
@printf("updated: %s\n", CONFIGFILE)
