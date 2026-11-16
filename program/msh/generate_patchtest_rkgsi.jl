# Generate the patch-test mesh with BenchmarkExample.PatchTest.generateMsh.
#
# generateMsh creates:
#   Ω, Γ¹..Γ⁴, and Γ (all cell edges via BenchmarkExample.addEdgeElements)

import BenchmarkExample: PatchTest
import BenchmarkExample.PatchTest: generateMsh

const n   = 10
const out = joinpath(@__DIR__, "patchtest_rkgsi.msh")

generateMsh(out; transfinite = n + 1, order = 1, quad = false)

println("written: ", out)
