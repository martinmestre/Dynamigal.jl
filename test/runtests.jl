using GalacticDynamics
using Test
using PythonCall


pyimport("sys")."path".append("")
pyimport("sys")."path".append("python")
accelerations_py = pyimport("accelerations")

@time begin
    include("acceleration/test_accelerations.jl")
end
@time begin
    include("orbit/test_orbits.jl")
end