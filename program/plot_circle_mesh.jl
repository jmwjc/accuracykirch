# ---------------------------------------------------------------------------
# plot_circle_mesh.jl
#
# Draw the full circular-plate mesh used by circle_rkgsi.jl.
#
# Run from the project root:
#
#     julia program/plot_circle_mesh.jl [ndiv]
#
# Output:
#     program/figures/circle_mesh.png
#     program/figures/circle_mesh_<ndiv>.png
#     manuscript/figs/circle_mesh.png
#     manuscript/figs/circle_mesh_<ndiv>.png
# ---------------------------------------------------------------------------

using CairoMakie
import Gmsh: gmsh

const R = 5.0
ndiv = parse(Int, get(ARGS, 1, "8"))
meshfile = joinpath(@__DIR__, "msh", "circle_full_ref_$(ndiv).msh")

function generate_full_disk(file::String, ndiv::Int)
    mkpath(dirname(file))
    gmsh.initialize()
    gmsh.model.add("Circle full")
    lc = R/ndiv
    gmsh.model.geo.addPoint(0.0, 0.0, 0.0, lc, 1)
    gmsh.model.geo.addPoint(R, 0.0, 0.0, lc, 2)
    gmsh.model.geo.addPoint(0.0, R, 0.0, lc, 3)
    gmsh.model.geo.addPoint(-R, 0.0, 0.0, lc, 4)
    gmsh.model.geo.addPoint(0.0, -R, 0.0, lc, 5)
    gmsh.model.geo.addCircleArc(2, 1, 3, 1)
    gmsh.model.geo.addCircleArc(3, 1, 4, 2)
    gmsh.model.geo.addCircleArc(4, 1, 5, 3)
    gmsh.model.geo.addCircleArc(5, 1, 2, 4)
    gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)
    gmsh.model.geo.addPlaneSurface([1], 1)
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.setAlgorithm(2, 1, 1)
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(1)
    gmsh.write(file)
    gmsh.finalize()
end

if !isfile(meshfile)
    generate_full_disk(meshfile, ndiv)
end

gmsh.initialize()
gmsh.open(meshfile)

nodeTags, coords, _ = gmsh.model.mesh.getNodes()
n = length(nodeTags)
xy = reshape(coords, 3, n)
points = [CairoMakie.Point2f(xy[1, i], xy[2, i]) for i in 1:n]
index = Dict(nodeTags[i] => i for i in 1:n)

elementTypes, _, elementNodeTags = gmsh.model.mesh.getElements(2)
faces = CairoMakie.TriangleFace{Int}[]
for (etype, ntags) in zip(elementTypes, elementNodeTags)
    etype == 2 || continue          # 3-node triangle
    for k in 1:3:length(ntags)
        push!(faces, CairoMakie.TriangleFace(
            index[ntags[k]], index[ntags[k + 1]], index[ntags[k + 2]]))
    end
end

gmsh.finalize()

fig = Figure(size = (700, 700))
ax = Axis(fig[1, 1], aspect = DataAspect())
line_strokewidth = ndiv == 4 ? 2.5 : ndiv == 8 ? 1.5 : ndiv == 16 ? 1.0 : 0.5
for f in faces
    poly!(ax, [points[f[1]], points[f[2]], points[f[3]]];
          color = :white, strokecolor = :black, strokewidth = line_strokewidth)
end

# Mesh nodes: filled black circles, smaller on finer meshes.
node_markersize = ndiv == 4 ? 32 : ndiv == 8 ? 20 : ndiv == 16 ? 12 : 8
node_mask = [hypot(p[1], p[2]) > 1e-10 for p in points]
scatter!(ax, points[node_mask];
         color = :black,
         marker = :circle,
         markersize = node_markersize,
         strokewidth = 0)
hidedecorations!(ax)
hidespines!(ax)

for dir in (joinpath(@__DIR__, "figures"),
            joinpath(@__DIR__, "..", "manuscript", "figs"))
    mkpath(dir)
    for name in ("circle_mesh.png", "circle_mesh_$(ndiv).png")
        save(joinpath(dir, name), fig; px_per_unit = 3)
        println("written: ", joinpath(dir, name))
    end
end
