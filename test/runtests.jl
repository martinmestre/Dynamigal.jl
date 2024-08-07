using GalacticDynamics
using Test
using PythonCall


pyimport("sys")."path".append("")
pyimport("sys")."path".append("python")
accelerations_py = pyimport("accelerations")
au = pyimport("astropy.units")
gd = pyimport("gala.dynamics")
gp = pyimport("gala.potential")
gu = pyimport("gala.units")

@time begin
    include("acceleration/test_accelerations.jl")
end
# @time begin
#     include("orbit/test_orbits.jl")
# end

