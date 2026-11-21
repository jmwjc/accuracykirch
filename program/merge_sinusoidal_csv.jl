# Merge p2/p3 XLSX results and the p4 CSV results into final CSV tables.
using XLSX

const XLSX_SOURCES = []
const CSV_BASE = "program/results/p4_support47.xlsx"
const OUTBASE = "program/sinusoidal_results"

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

function read_xlsx_table(path, sheet)
    t = XLSX.readtable(path, sheet)
    return t.column_labels, [Any[x for x in col] for col in t.data]
end

function append_cols!(dst, src)
    for i in eachindex(dst)
        append!(dst[i], src[i])
    end
end

function csv_cell(x)
    s = ismissing(x) ? "" : string(x)
    if occursin(',', s) || occursin('"', s) || occursin('\n', s)
        return "\"" * replace(s, "\"" => "\"\"") * "\""
    end
    return s
end

function write_csv(file, labels, cols)
    open(file, "w") do io
        println(io, join(string.(labels), ","))
        n = length(cols[1])
        for i in 1:n
            println(io, join([csv_cell(cols[j][i]) for j in eachindex(cols)], ","))
        end
    end
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

function merge_sheet(sheet)
    labels = nothing
    cols = nothing
    for f in XLSX_SOURCES
        isfile(f) || continue
        l, c = read_xlsx_table(f, sheet)
        labels === nothing && (labels = l; cols = [Any[] for _ in l])
        append_cols!(cols, c)
    end
    csvfile = sheet == "raw" ? CSV_BASE * ".csv" :
              sheet == "rates" ? CSV_BASE * "_rates.csv" :
                                 CSV_BASE * "_config.csv"
    if isfile(csvfile)
        l, c = read_csv_table(csvfile)
        labels === nothing && (labels = l; cols = [Any[] for _ in l])
        append_cols!(cols, c)
    end
    return labels, cols
end

lraw, craw = merge_sheet("raw")
iv = findfirst(==(:valid), lraw)
if iv !== nothing
    for j in eachindex(craw[1])
        if craw[iv][j] == "no"
            for name in (:L2, :H1, :H2, :He)
                k = findfirst(==(name), lraw)
                k === nothing || (craw[k][j] = "")
            end
        end
    end
end
sort_rows!(lraw, craw, (:p, :q, :ndiv))

lrate, crate = merge_sheet("rates")
sort_rows!(lrate, crate, (:p, :q, :norm), ["L2", "H1", "H2"])

lconf, cconf = merge_sheet("config")
dedup_config!(lconf, cconf)

write_csv(OUTBASE * ".csv", lraw, craw)
write_csv(OUTBASE * "_rates.csv", lrate, crate)
write_csv(OUTBASE * "_config.csv", lconf, cconf)
println("written: ", OUTBASE * ".csv")
println("written: ", OUTBASE * "_rates.csv")
println("written: ", OUTBASE * "_config.csv")
